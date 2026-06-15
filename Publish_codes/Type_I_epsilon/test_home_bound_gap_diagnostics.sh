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
B_VALUES="${B_VALUES:-0.025 0.1 0.2 1 10}"
QUERY_INDEXES="${QUERY_INDEXES:-0 10 100}"
OUTPUT_ROOT="${OUTPUT_ROOT:-${ROOT_DIR}/test_Result/home_bound_gap_diagnostics}"
BIN="${BIN:-${ROOT_DIR}/build/release/bound_diagnostics}"

if [[ ! -x "${BIN}" ]]; then
	echo "Error: diagnostic binary not found: ${BIN}" >&2
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

summarize_diag() {
	local b="$1"
	local query_index="$2"
	local csv="$3"
	awk -F, -v b="${b}" -v q="${query_index}" '
	NR==1 {next}
	{
		nodes++;
		dx=$7+0;
		du=$8+0;
		z=$9+0;
		r=$10+0;
		k_gap=$13+0;
		k_rel=$14+0;
		valid=$15+0;
		a_gap=$18+0;
		a_rel=$19+0;
		a_contains=$24+0;
		a_abs_better=$25+0;
		a_rel_better=$26+0;
		zr=(r>0 ? z/r : -1);
		dd=(dx>0 ? du/dx : -1);
		if(!valid) invalid++;
		if(valid && !a_contains) anchor_miss++;
		if(a_abs_better) abs_better++;
		if(a_rel_better) {
			rel_better++;
			rel_dx+=dx;
			rel_du+=du;
			rel_dd+=dd;
			rel_zr+=zr;
			rel_k_gap+=k_gap;
			rel_a_gap+=a_gap;
			if(zr<=1) rel_near++;
			else if(zr<=3) rel_mid++;
			else rel_far++;
		} else {
			karl_better++;
			kar_dx+=dx;
			kar_du+=du;
			kar_dd+=dd;
			kar_zr+=zr;
			kar_k_gap+=k_gap;
			kar_a_gap+=a_gap;
			if(zr<=1) kar_near++;
			else if(zr<=3) kar_mid++;
			else kar_far++;
		}
	}
	END {
		printf "%s,%s,%d,%d,%.6f,%d,%.6f,%d,%d,", b,q,nodes,abs_better,(nodes?abs_better/nodes:0),rel_better,(nodes?rel_better/nodes:0),invalid,anchor_miss;
		printf "%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,", (rel_better?rel_dx/rel_better:0), (rel_better?rel_du/rel_better:0), (rel_better?rel_dd/rel_better:0), (rel_better?rel_zr/rel_better:0), (rel_better?rel_k_gap/rel_better:0), (rel_better?rel_a_gap/rel_better:0);
		printf "%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,", (karl_better?kar_dx/karl_better:0), (karl_better?kar_du/karl_better:0), (karl_better?kar_dd/karl_better:0), (karl_better?kar_zr/karl_better:0), (karl_better?kar_k_gap/karl_better:0), (karl_better?kar_a_gap/karl_better:0);
		printf "%d,%d,%d,%d,%d,%d\n", rel_near,rel_mid,rel_far,kar_near,kar_mid,kar_far;
	}' "${csv}"
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

echo "Diagnostic binary: ${BIN}"
echo "Query subset: ${QUERY_FILE}"
echo "Data subset: ${DATA_FILE}"
echo "B_VALUES=${B_VALUES}"
echo "QUERY_INDEXES=${QUERY_INDEXES}"
echo "Run dir: ${RUN_DIR}"
echo

printf '%s\n' \
	"b,query_index,nodes,anchor_abs_better,anchor_abs_better_frac,anchor_rel_better,anchor_rel_better_frac,anchor_invalid,anchor_contains_fail,anchor_better_avg_delta_x,anchor_better_avg_delta_u,anchor_better_avg_du_over_dx,anchor_better_avg_z_over_radius,anchor_better_avg_karl_gap,anchor_better_avg_anchor_gap,karl_better_avg_delta_x,karl_better_avg_delta_u,karl_better_avg_du_over_dx,karl_better_avg_z_over_radius,karl_better_avg_karl_gap,karl_better_avg_anchor_gap,anchor_better_near,anchor_better_mid,anchor_better_far,karl_better_near,karl_better_mid,karl_better_far" \
	> "${SUMMARY_CSV}"

for B in ${B_VALUES}; do
	B_TAG="${B//./p}"
	for Q_INDEX in ${QUERY_INDEXES}; do
		OUT_CSV="${RUN_DIR}/diag_b${B_TAG}_q${Q_INDEX}.csv"
		echo "[RUN] b=${B} query_index=${Q_INDEX}"
		"${BIN}" "${QUERY_FILE}" "${DATA_FILE}" "${OUT_CSV}" "${LEAF_CAPACITY}" "${REL_ERROR}" "${B}" "${Q_INDEX}"
		summarize_diag "${B}" "${Q_INDEX}" "${OUT_CSV}" >> "${SUMMARY_CSV}"
	done
done

echo
echo "Summary CSV: ${SUMMARY_CSV}"
cat "${SUMMARY_CSV}"
