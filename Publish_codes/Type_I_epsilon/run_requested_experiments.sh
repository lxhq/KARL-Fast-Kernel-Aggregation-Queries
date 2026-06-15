#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="${SCRIPT_DIR}"

B_VALUE="${B_VALUE:-0.1}"
REL_ERROR="${REL_ERROR:-0.2}"
OUTPUT_ROOT="${OUTPUT_ROOT:-${ROOT_DIR}/test_Result/requested_experiments}"

if [[ -x "${ROOT_DIR}/build/release/main_release" ]]; then
    BIN="${BIN:-${ROOT_DIR}/build/release/main_release}"
elif [[ -x "${ROOT_DIR}/main" ]]; then
    BIN="${BIN:-${ROOT_DIR}/main}"
elif [[ -x "${ROOT_DIR}/build/debug/main_debug" ]]; then
    BIN="${BIN:-${ROOT_DIR}/build/debug/main_debug}"
else
    echo "Error: no runnable binary found." >&2
    echo "Expected one of: ./build/release/main_release, ./main, ./build/debug/main_debug" >&2
    exit 1
fi

DATA_ROOT="/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation"

SUSY_Q="${DATA_ROOT}/susy/SUSY_qSet.data"
SUSY_X="${DATA_ROOT}/susy/SUSY_X.data"
HOME_Q="${DATA_ROOT}/home/HT_Sensor_dataset_qSet.data"
HOME_X="${DATA_ROOT}/home/HT_Sensor_dataset_X.data"
MINI_Q="${DATA_ROOT}/miniboone/MiniBooNE_qSet.data"
MINI_X="${DATA_ROOT}/miniboone/MiniBooNE_X.data"

require_file() {
    local path="$1"
    if [[ ! -f "${path}" ]]; then
        echo "Error: required file not found: ${path}" >&2
        exit 1
    fi
}

query_rows_from_matrix() {
    local matrix_path="$1"
    awk 'NR == 1 {print $1; exit}' "${matrix_path}"
}

mkdir -p "${OUTPUT_ROOT}"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${OUTPUT_ROOT}/run_${RUN_ID}"
mkdir -p "${RUN_DIR}"

SUMMARY_CSV="${RUN_DIR}/summary.csv"
SUMMARY_TXT="${RUN_DIR}/summary.txt"

printf 'case_id,dataset,method,leaf,b,rel_error,query_count,total_time_s,qps,status,query_file,dataset_file,result_file,log_file\n' > "${SUMMARY_CSV}"
printf '%-30s %-10s %-6s %-6s %-6s %-10s %-10s %-12s %-12s %-10s\n' \
    "case_id" "dataset" "method" "leaf" "b" "rel_error" "q_count" "time_s" "qps" "status" > "${SUMMARY_TXT}"

echo "Binary: ${BIN}"
echo "b: ${B_VALUE}"
echo "rel_error: ${REL_ERROR}"
echo "Run dir: ${RUN_DIR}"
echo

for f in "${SUSY_Q}" "${SUSY_X}" "${HOME_Q}" "${HOME_X}" "${MINI_Q}" "${MINI_X}"; do
    require_file "${f}"
done

run_case() {
    local case_id="$1"
    local dataset="$2"
    local query_file="$3"
    local dataset_file="$4"
    local method="$5"
    local leaf="$6"

    local result_file="${RUN_DIR}/${case_id}.txt"
    local log_file="${RUN_DIR}/${case_id}.log"
    local query_count
    local total_time
    local qps
    local status="ok"

    query_count="$(query_rows_from_matrix "${query_file}")"

    echo "[RUN] ${case_id}"
    echo "      ${BIN} ${query_file} ${dataset_file} ${result_file} ${method} ${leaf} ${REL_ERROR} ${B_VALUE}"

    set +e
    "${BIN}" "${query_file}" "${dataset_file}" "${result_file}" "${method}" "${leaf}" "${REL_ERROR}" "${B_VALUE}" \
        2>&1 | tee "${log_file}"
    local rc=${PIPESTATUS[0]}
    set -e

    if (( rc != 0 )); then
        status="failed(${rc})"
    fi

    total_time="$(awk '/^Total time:/ {print $3}' "${log_file}" | tail -n 1)"
    qps="$(awk -v m="${method}" '$1=="Method" {id=$2; gsub(":", "", id); if (id==m) print $3}' "${log_file}" | tail -n 1)"

    if [[ -z "${total_time}" ]]; then
        total_time="NA"
    fi
    if [[ -z "${qps}" ]]; then
        qps="NA"
    fi

    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "${case_id}" "${dataset}" "${method}" "${leaf}" "${B_VALUE}" "${REL_ERROR}" "${query_count}" \
        "${total_time}" "${qps}" "${status}" "${query_file}" "${dataset_file}" "${result_file}" "${log_file}" \
        >> "${SUMMARY_CSV}"

    printf '%-30s %-10s %-6s %-6s %-6s %-10s %-10s %-12s %-12s %-10s\n' \
        "${case_id}" "${dataset}" "${method}" "${leaf}" "${B_VALUE}" "${REL_ERROR}" "${query_count}" \
        "${total_time}" "${qps}" "${status}" >> "${SUMMARY_TXT}"

    echo
}

# Requested method 0 runs (with qSet + X paths from launch.json)
run_case "susy_m0_leaf10_b0p1" "susy" "${SUSY_Q}" "${SUSY_X}" "0" "10"
run_case "home_m0_leaf10_b0p1" "home" "${HOME_Q}" "${HOME_X}" "0" "10"
run_case "miniboone_m0_leaf10_b0p1" "miniboone" "${MINI_Q}" "${MINI_X}" "0" "10"

# Additional requested runs
run_case "home_m3_leaf80_b0p1" "home" "${HOME_Q}" "${HOME_X}" "3" "80"
run_case "miniboone_m12_leaf10_b0p1" "miniboone" "${MINI_Q}" "${MINI_X}" "12" "10"
run_case "susy_m3_leaf10_b0p1" "susy" "${SUSY_Q}" "${SUSY_X}" "3" "10"

echo "Finished."
echo "Summary CSV: ${SUMMARY_CSV}"
echo "Summary TXT: ${SUMMARY_TXT}"
