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
OUTPUT_ROOT="${OUTPUT_ROOT:-${ROOT_DIR}/test_Result/home_anchor_bandwidth_profile}"
B_VALUES="${B_VALUES:-0.025 0.05 0.1 0.2 0.5 1 2 5 10}"

if [[ -x "${ROOT_DIR}/build/release/main_release" ]]; then
	BIN="${BIN:-${ROOT_DIR}/build/release/main_release}"
elif [[ -x "${ROOT_DIR}/main" ]]; then
	BIN="${BIN:-${ROOT_DIR}/main}"
else
	echo "Error: no runnable KARL binary found." >&2
	exit 1
fi

require_file() {
	local path="$1"
	if [[ ! -f "${path}" ]]; then
		echo "Error: required file not found: ${path}" >&2
		exit 1
	fi
}

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

extract_total_time() {
	local log_file="$1"
	awk '/^Total time:/ {print $3}' "${log_file}" | tail -n 1
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
		abs_err_sum+=diff;
		if(diff>max_abs_err) max_abs_err=diff;
		denom=( (a<0?-a:a) + (r<0?-r:r) + 1e-12 );
		sym_rel=(2.0*diff)/denom;
		sym_rel_sum+=sym_rel;
		if(sym_rel>max_sym_rel) max_sym_rel=sym_rel;
		valid_pairs++;
	}
	END {
		if(valid_pairs==0) {
			print "nan,nan,nan,nan,0," nan_pairs;
		} else {
			printf "%.10g,%.10g,%.10g,%.10g,%d,%d\n",
				abs_err_sum/valid_pairs, max_abs_err,
				sym_rel_sum/valid_pairs, max_sym_rel,
				valid_pairs, nan_pairs;
		}
	}' "${ref_file}" "${approx_file}"
}

require_file "${QUERY_SOURCE}"
require_file "${DATA_SOURCE}"

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
echo "b values: ${B_VALUES}"
echo "Run dir: ${RUN_DIR}"
echo

printf '%s\n' \
	"b,m0_qps,m3_qps,m13_qps,m14_qps,m13_over_m3,m14_over_m3,m3_speedup_vs_m0,m13_speedup_vs_m0,m14_speedup_vs_m0,m3_avg_nodes,m13_avg_nodes,m14_avg_nodes,m3_p95_nodes,m13_p95_nodes,m14_p95_nodes,m3_exact_ratio,m13_exact_ratio,m14_exact_ratio,m3_avg_bound_calls,m13_avg_bound_calls,m14_avg_bound_calls,m3_time_bound_sec,m13_time_bound_sec,m14_time_bound_sec,m3_mean_sym_rel,m13_mean_sym_rel,m14_mean_sym_rel,m3_valid_pairs,m13_valid_pairs,m14_valid_pairs" \
	> "${SUMMARY_CSV}"

for B in ${B_VALUES}; do
	B_TAG="${B//./p}"
	B_TAG="${B_TAG//-/m}"
	OUT_M0="${RUN_DIR}/home_b${B_TAG}_m0.txt"
	OUT_M3="${RUN_DIR}/home_b${B_TAG}_m3.txt"
	OUT_M13="${RUN_DIR}/home_b${B_TAG}_m13.txt"
	OUT_M14="${RUN_DIR}/home_b${B_TAG}_m14.txt"
	LOG_M0="${RUN_DIR}/home_b${B_TAG}_m0.log"
	LOG_M3="${RUN_DIR}/home_b${B_TAG}_m3.log"
	LOG_M13="${RUN_DIR}/home_b${B_TAG}_m13.log"
	LOG_M14="${RUN_DIR}/home_b${B_TAG}_m14.log"

	echo "[RUN] b=${B} method=0"
	"${BIN}" "${QUERY_FILE}" "${DATA_FILE}" "${OUT_M0}" 0 "${LEAF_CAPACITY}" "${REL_ERROR}" "${B}" | tee "${LOG_M0}"

	echo "[RUN] b=${B} method=3"
	"${BIN}" "${QUERY_FILE}" "${DATA_FILE}" "${OUT_M3}" 3 "${LEAF_CAPACITY}" "${REL_ERROR}" "${B}" | tee "${LOG_M3}"

	echo "[RUN] b=${B} method=13"
	"${BIN}" "${QUERY_FILE}" "${DATA_FILE}" "${OUT_M13}" 13 "${LEAF_CAPACITY}" "${REL_ERROR}" "${B}" | tee "${LOG_M13}"

	echo "[RUN] b=${B} method=14"
	"${BIN}" "${QUERY_FILE}" "${DATA_FILE}" "${OUT_M14}" 14 "${LEAF_CAPACITY}" "${REL_ERROR}" "${B}" | tee "${LOG_M14}"

	M0_QPS="$(extract_qps 0 "${LOG_M0}")"
	M3_QPS="$(extract_qps 3 "${LOG_M3}")"
	M13_QPS="$(extract_qps 13 "${LOG_M13}")"
	M14_QPS="$(extract_qps 14 "${LOG_M14}")"
	M13_OVER_M3="$(awk -v a="${M13_QPS}" -v b="${M3_QPS}" 'BEGIN{if(b==0) print "nan"; else printf "%.6f", a/b}')"
	M14_OVER_M3="$(awk -v a="${M14_QPS}" -v b="${M3_QPS}" 'BEGIN{if(b==0) print "nan"; else printf "%.6f", a/b}')"
	M3_OVER_M0="$(awk -v a="${M3_QPS}" -v b="${M0_QPS}" 'BEGIN{if(b==0) print "nan"; else printf "%.6f", a/b}')"
	M13_OVER_M0="$(awk -v a="${M13_QPS}" -v b="${M0_QPS}" 'BEGIN{if(b==0) print "nan"; else printf "%.6f", a/b}')"
	M14_OVER_M0="$(awk -v a="${M14_QPS}" -v b="${M0_QPS}" 'BEGIN{if(b==0) print "nan"; else printf "%.6f", a/b}')"

	M3_AVG_NODES="$(extract_metric "avg_nodes_processed_per_query" "${LOG_M3}")"
	M13_AVG_NODES="$(extract_metric "avg_nodes_processed_per_query" "${LOG_M13}")"
	M14_AVG_NODES="$(extract_metric "avg_nodes_processed_per_query" "${LOG_M14}")"
	M3_P95_NODES="$(extract_metric "p95_nodes_processed_per_query" "${LOG_M3}")"
	M13_P95_NODES="$(extract_metric "p95_nodes_processed_per_query" "${LOG_M13}")"
	M14_P95_NODES="$(extract_metric "p95_nodes_processed_per_query" "${LOG_M14}")"
	M3_EXACT_RATIO="$(extract_metric "exact_point_eval_ratio_vs_full_scan" "${LOG_M3}")"
	M13_EXACT_RATIO="$(extract_metric "exact_point_eval_ratio_vs_full_scan" "${LOG_M13}")"
	M14_EXACT_RATIO="$(extract_metric "exact_point_eval_ratio_vs_full_scan" "${LOG_M14}")"
	M3_BOUND_CALLS="$(extract_metric "avg_bound_calls_per_query" "${LOG_M3}")"
	M13_BOUND_CALLS="$(extract_metric "avg_bound_calls_per_query" "${LOG_M13}")"
	M14_BOUND_CALLS="$(extract_metric "avg_bound_calls_per_query" "${LOG_M14}")"
	M3_TIME_BOUND="$(extract_metric "time_bound_sec" "${LOG_M3}")"
	M13_TIME_BOUND="$(extract_metric "time_bound_sec" "${LOG_M13}")"
	M14_TIME_BOUND="$(extract_metric "time_bound_sec" "${LOG_M14}")"

	M3_ERRORS="$(compare_outputs "${OUT_M0}" "${OUT_M3}")"
	M13_ERRORS="$(compare_outputs "${OUT_M0}" "${OUT_M13}")"
	M14_ERRORS="$(compare_outputs "${OUT_M0}" "${OUT_M14}")"
	M3_MEAN_SYM="$(awk -F, '{print $3}' <<< "${M3_ERRORS}")"
	M13_MEAN_SYM="$(awk -F, '{print $3}' <<< "${M13_ERRORS}")"
	M14_MEAN_SYM="$(awk -F, '{print $3}' <<< "${M14_ERRORS}")"
	M3_VALID="$(awk -F, '{print $5}' <<< "${M3_ERRORS}")"
	M13_VALID="$(awk -F, '{print $5}' <<< "${M13_ERRORS}")"
	M14_VALID="$(awk -F, '{print $5}' <<< "${M14_ERRORS}")"

	printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
		"${B}" "${M0_QPS}" "${M3_QPS}" "${M13_QPS}" "${M14_QPS}" \
		"${M13_OVER_M3}" "${M14_OVER_M3}" "${M3_OVER_M0}" "${M13_OVER_M0}" "${M14_OVER_M0}" \
		"${M3_AVG_NODES}" "${M13_AVG_NODES}" "${M14_AVG_NODES}" \
		"${M3_P95_NODES}" "${M13_P95_NODES}" "${M14_P95_NODES}" \
		"${M3_EXACT_RATIO}" "${M13_EXACT_RATIO}" "${M14_EXACT_RATIO}" \
		"${M3_BOUND_CALLS}" "${M13_BOUND_CALLS}" "${M14_BOUND_CALLS}" \
		"${M3_TIME_BOUND}" "${M13_TIME_BOUND}" "${M14_TIME_BOUND}" \
		"${M3_MEAN_SYM}" "${M13_MEAN_SYM}" "${M14_MEAN_SYM}" \
		"${M3_VALID}" "${M13_VALID}" "${M14_VALID}" >> "${SUMMARY_CSV}"

	echo "[DONE] b=${B} m13/m3=${M13_OVER_M3} m14/m3=${M14_OVER_M3} m3_nodes=${M3_AVG_NODES} m13_nodes=${M13_AVG_NODES} m14_nodes=${M14_AVG_NODES}"
	echo
done

echo "Summary CSV: ${SUMMARY_CSV}"
cat "${SUMMARY_CSV}"
