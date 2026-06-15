#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/home"

QUERY_SOURCE="${QUERY_SOURCE:-${DATA_ROOT}/HT_Sensor_dataset_sampleSet.data}"
DATA_SOURCE="${DATA_SOURCE:-${DATA_ROOT}/HT_Sensor_dataset_X.data}"
Q_NUM="${Q_NUM:-1000}"
X_NUM="${X_NUM:-100000}"
DIM="${DIM:-10}"
REL_ERROR="${REL_ERROR:-0.2}"
LEAF_CAPACITY="${LEAF_CAPACITY:-80}"
B_VALUES="${B_VALUES:-0.025 0.05 0.1 0.2 0.5 1 2 5 10}"
OUTPUT_ROOT="${OUTPUT_ROOT:-${ROOT_DIR}/test_Result/home_adaptive_anchor_profile}"
BIN="${BIN:-${ROOT_DIR}/build/release/main_release}"

if [[ ! -x "${BIN}" ]]; then
	echo "Error: binary not found: ${BIN}" >&2
	exit 1
fi

make_subset() {
	local src="$1"
	local dst="$2"
	local rows="$3"
	local dim="$4"

	{
		printf '%s %s\n' "${rows}" "${dim}"
		sed -n "2,$((rows+1))p" "${src}"
	} > "${dst}"
}

extract_qps() {
	local method="$1"
	local log_file="$2"
	awk -v method="${method}" '$1 == "Method" {
		id=$2;
		gsub(":", "", id);
		if(id == method) print $3;
	}' "${log_file}" | tail -n 1
}

extract_metric() {
	local metric="$1"
	local log_file="$2"
	awk -v metric="${metric}" '
	{
		n=split($0, parts, /[, ]+/);
		for(i=1; i<=n; i++) {
			if(index(parts[i], metric "=") == 1) {
				split(parts[i], kv, "=");
				print kv[2];
				exit;
			}
		}
	}' "${log_file}"
}

compare_outputs() {
	local ref_file="$1"
	local approx_file="$2"
	awk '
	NR==FNR {
		ref[NR]=$1;
		next;
	}
	{
		r=ref[NR];
		a=$1;
		if(tolower(r) ~ /nan/ || tolower(a) ~ /nan/) {
			nan_pairs++;
			next;
		}
		diff=a-r;
		if(diff<0) diff=-diff;
		denom=( (a<0?-a:a) + (r<0?-r:r) + 1e-12 );
		sym_rel=(2.0*diff)/denom;
		sym_rel_sum+=sym_rel;
		valid_pairs++;
	}
	END {
		if(valid_pairs==0) {
			print "nan,0," nan_pairs;
		} else {
			printf "%.10g,%d,%d\n", sym_rel_sum/valid_pairs, valid_pairs, nan_pairs;
		}
	}' "${ref_file}" "${approx_file}"
}

mkdir -p "${OUTPUT_ROOT}"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${OUTPUT_ROOT}/run_${RUN_ID}"
mkdir -p "${RUN_DIR}"

QUERY_FILE="${RUN_DIR}/home_q${Q_NUM}.data"
DATA_FILE="${RUN_DIR}/home_x${X_NUM}.data"
SUMMARY_CSV="${RUN_DIR}/summary.csv"
MASTER_LOG="${RUN_DIR}/run.log"

exec > >(tee -a "${MASTER_LOG}") 2>&1

make_subset "${QUERY_SOURCE}" "${QUERY_FILE}" "${Q_NUM}" "${DIM}"
make_subset "${DATA_SOURCE}" "${DATA_FILE}" "${X_NUM}" "${DIM}"

echo "Binary: ${BIN}"
echo "Query subset: ${QUERY_FILE}"
echo "Data subset: ${DATA_FILE}"
echo "Q_NUM=${Q_NUM}, X_NUM=${X_NUM}, DIM=${DIM}, REL_ERROR=${REL_ERROR}, LEAF_CAPACITY=${LEAF_CAPACITY}"
echo "B_VALUES=${B_VALUES}"
echo "Run dir: ${RUN_DIR}"
echo

printf '%s\n' \
	"b,m0_qps,m3_qps,m14_qps,m15_qps,m16_qps,m14_over_m3,m15_over_m3,m16_over_m3,m3_nodes,m14_nodes,m15_nodes,m16_nodes,m3_exact_ratio,m14_exact_ratio,m15_exact_ratio,m16_exact_ratio,m15_mean_sym_rel,m16_mean_sym_rel,m15_valid_pairs,m16_valid_pairs" \
	> "${SUMMARY_CSV}"

for B in ${B_VALUES}; do
	B_TAG="${B//./p}"
	declare -A OUT
	declare -A LOG
	declare -A QPS
	declare -A NODES
	declare -A EXACT_RATIO

	for METHOD in 0 3 14 15 16; do
		OUT[${METHOD}]="${RUN_DIR}/home_b${B_TAG}_m${METHOD}.txt"
		LOG[${METHOD}]="${RUN_DIR}/home_b${B_TAG}_m${METHOD}.log"
		echo "[RUN] b=${B} method=${METHOD}"
		"${BIN}" "${QUERY_FILE}" "${DATA_FILE}" "${OUT[${METHOD}]}" "${METHOD}" "${LEAF_CAPACITY}" "${REL_ERROR}" "${B}" | tee "${LOG[${METHOD}]}"
		QPS[${METHOD}]="$(extract_qps "${METHOD}" "${LOG[${METHOD}]}")"
		if [[ "${METHOD}" != "0" ]]; then
			NODES[${METHOD}]="$(extract_metric "avg_nodes_processed_per_query" "${LOG[${METHOD}]}")"
			EXACT_RATIO[${METHOD}]="$(extract_metric "exact_point_eval_ratio_vs_full_scan" "${LOG[${METHOD}]}")"
		fi
	done

	M14_OVER_M3="$(awk -v a="${QPS[14]}" -v b="${QPS[3]}" 'BEGIN{if(b==0) print "nan"; else printf "%.6f", a/b}')"
	M15_OVER_M3="$(awk -v a="${QPS[15]}" -v b="${QPS[3]}" 'BEGIN{if(b==0) print "nan"; else printf "%.6f", a/b}')"
	M16_OVER_M3="$(awk -v a="${QPS[16]}" -v b="${QPS[3]}" 'BEGIN{if(b==0) print "nan"; else printf "%.6f", a/b}')"
	M15_ERR="$(compare_outputs "${OUT[0]}" "${OUT[15]}")"
	M16_ERR="$(compare_outputs "${OUT[0]}" "${OUT[16]}")"
	M15_MEAN="$(awk -F, '{print $1}' <<< "${M15_ERR}")"
	M16_MEAN="$(awk -F, '{print $1}' <<< "${M16_ERR}")"
	M15_VALID="$(awk -F, '{print $2}' <<< "${M15_ERR}")"
	M16_VALID="$(awk -F, '{print $2}' <<< "${M16_ERR}")"

	printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
		"${B}" "${QPS[0]}" "${QPS[3]}" "${QPS[14]}" "${QPS[15]}" "${QPS[16]}" \
		"${M14_OVER_M3}" "${M15_OVER_M3}" "${M16_OVER_M3}" \
		"${NODES[3]}" "${NODES[14]}" "${NODES[15]}" "${NODES[16]}" \
		"${EXACT_RATIO[3]}" "${EXACT_RATIO[14]}" "${EXACT_RATIO[15]}" "${EXACT_RATIO[16]}" \
		"${M15_MEAN}" "${M16_MEAN}" "${M15_VALID}" "${M16_VALID}" >> "${SUMMARY_CSV}"

	echo "[DONE] b=${B} m14/m3=${M14_OVER_M3} m15/m3=${M15_OVER_M3} m16/m3=${M16_OVER_M3}"
	echo
done

echo "Summary CSV: ${SUMMARY_CSV}"
cat "${SUMMARY_CSV}"
