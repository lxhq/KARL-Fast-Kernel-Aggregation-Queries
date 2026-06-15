#!/usr/bin/env bash
set -euo pipefail

# Set your test inputs here (can be overridden by environment variables).
QUERY_FILE="${QUERY_FILE:-/tmp/mini_q300.data}"
DATASET_FILE="${DATASET_FILE:-/tmp/mini_x30000.data}"
REL_ERROR="${REL_ERROR:-0.2}"
LEAF_CAPACITY="${LEAF_CAPACITY:-40}"
OUTPUT_DIR="${OUTPUT_DIR:-./test_Result/b_sweep}"

if [[ -n "${B_VALUES_ENV:-}" ]]; then
    read -r -a B_VALUES <<< "$B_VALUES_ENV"
else
    B_VALUES=(0.5 1 2 5 10 20)
fi

if [[ ! -f "$QUERY_FILE" ]]; then
    echo "Error: query file not found: $QUERY_FILE" >&2
    exit 1
fi

if [[ ! -f "$DATASET_FILE" ]]; then
    echo "Error: dataset file not found: $DATASET_FILE" >&2
    exit 1
fi

if [[ -x "./build/release/main_release" ]]; then
    BIN="./build/release/main_release"
elif [[ -x "./main" ]]; then
    BIN="./main"
elif [[ -x "./build/debug/main_debug" ]]; then
    BIN="./build/debug/main_debug"
else
    echo "Error: no runnable binary found." >&2
    echo "Expected one of: ./build/release/main_release, ./main, ./build/debug/main_debug" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${OUTPUT_DIR}/run_${TIMESTAMP}"
mkdir -p "$RUN_DIR"
SUMMARY_CSV="${RUN_DIR}/summary.csv"
MASTER_LOG="${RUN_DIR}/run.log"

exec > >(tee -a "$MASTER_LOG") 2>&1

echo "Binary: $BIN"
echo "Query file: $QUERY_FILE"
echo "Dataset file: $DATASET_FILE"
echo "Relative error: $REL_ERROR"
echo "Leaf capacity: $LEAF_CAPACITY"
echo "b values: ${B_VALUES[*]}"
echo "Run dir: $RUN_DIR"
echo

echo "b,m0_qps,m3_qps,m3_over_m0,mae,max_abs_err,mean_sym_rel_err,max_sym_rel_err,valid_pairs,nan_pairs" > "$SUMMARY_CSV"

for B in "${B_VALUES[@]}"; do
    B_TAG="${B//./p}"
    OUT_M0="${RUN_DIR}/method0_b_${B_TAG}.txt"
    OUT_M3="${RUN_DIR}/method3_b_${B_TAG}.txt"
    LOG_M0="${RUN_DIR}/method0_b_${B_TAG}.log"
    LOG_M3="${RUN_DIR}/method3_b_${B_TAG}.log"

    echo "[RUN] b=${B} method 0"
    "$BIN" "$QUERY_FILE" "$DATASET_FILE" "$OUT_M0" 0 "$LEAF_CAPACITY" "$REL_ERROR" "$B" | tee "$LOG_M0"

    echo "[RUN] b=${B} method 3"
    "$BIN" "$QUERY_FILE" "$DATASET_FILE" "$OUT_M3" 3 "$LEAF_CAPACITY" "$REL_ERROR" "$B" | tee "$LOG_M3"

    M0_QPS="$(awk '/^Method 0:/{print $3}' "$LOG_M0" | tail -n 1)"
    M3_QPS="$(awk '/^Method 3:/{print $3}' "$LOG_M3" | tail -n 1)"

    read -r MAE MAX_ABS_ERR MEAN_SYM_REL MAX_SYM_REL VALID_PAIRS NAN_PAIRS <<EOF
$(awk '
NR==FNR { ref[NR]=$1; next }
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
        print "nan nan nan nan 0 " nan_pairs;
    } else {
        printf "%.10g %.10g %.10g %.10g %d %d\n", abs_err_sum/valid_pairs, max_abs_err, sym_rel_sum/valid_pairs, max_sym_rel, valid_pairs, nan_pairs;
    }
}' "$OUT_M0" "$OUT_M3")
EOF

    M3_OVER_M0="$(awk -v a="$M3_QPS" -v b="$M0_QPS" 'BEGIN{if(b==0) print "nan"; else printf "%.6f", a/b}')"

    echo "${B},${M0_QPS},${M3_QPS},${M3_OVER_M0},${MAE},${MAX_ABS_ERR},${MEAN_SYM_REL},${MAX_SYM_REL},${VALID_PAIRS},${NAN_PAIRS}" >> "$SUMMARY_CSV"
    echo "[DONE] b=${B} valid_pairs=${VALID_PAIRS} nan_pairs=${NAN_PAIRS} m3/m0=${M3_OVER_M0} mae=${MAE} mean_sym_rel=${MEAN_SYM_REL}"
    echo
done

echo "Summary CSV: $SUMMARY_CSV"
cat "$SUMMARY_CSV"
