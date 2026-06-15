#!/usr/bin/env bash
set -euo pipefail

BIN="./build/release/main_release"
if [[ ! -x "$BIN" ]]; then
    echo "Error: binary not found: $BIN" >&2
    exit 1
fi

LEAF_CAPACITY="${LEAF_CAPACITY:-40}"
REL_ERROR="${REL_ERROR:-0.2}"
B_VALUE="${B_VALUE:-1}"
SCALES="${SCALES:-0.25 0.5 1 2 4}"
OUT_ROOT="${OUT_ROOT:-./test_Result/scale_effect}"

HT_Q="/tmp/karl_exp_ht_q300.data"
HT_X="/tmp/karl_exp_ht_x30000.data"
MINI_Q="/tmp/karl_exp_miniboone_q300.data"
MINI_X="/tmp/karl_exp_miniboone_x30000.data"

for f in "$HT_Q" "$HT_X" "$MINI_Q" "$MINI_X"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: required subset file missing: $f" >&2
        exit 1
    fi
done

mkdir -p "$OUT_ROOT"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${OUT_ROOT}/run_${RUN_ID}"
mkdir -p "$RUN_DIR"
SUMMARY_CSV="${RUN_DIR}/summary.csv"

echo "dataset,scale,m0_qps,m3_qps,m3_over_m0,raw_mean_abs_q,raw_mean_abs_x" > "$SUMMARY_CSV"

make_scaled_file() {
    local in_file="$1"
    local out_file="$2"
    local s="$3"
    awk -v s="$s" 'NR==1{print;next}{for(i=1;i<=NF;i++) $i=$i*s; print}' "$in_file" > "$out_file"
}

mean_abs_file() {
    local f="$1"
    awk 'NR==1{next}{for(i=1;i<=NF;i++){v=$i; if(v<0)v=-v; sum+=v; n++}} END{if(n==0)print "nan"; else printf "%.6f",sum/n}' "$f"
}

run_one() {
    local dataset="$1"
    local q_file="$2"
    local x_file="$3"
    local scale="$4"

    local tag_scale="${scale//./p}"
    local out0="${RUN_DIR}/${dataset}_s${tag_scale}_m0.txt"
    local out3="${RUN_DIR}/${dataset}_s${tag_scale}_m3.txt"
    local log0="${RUN_DIR}/${dataset}_s${tag_scale}_m0.log"
    local log3="${RUN_DIR}/${dataset}_s${tag_scale}_m3.log"

    "$BIN" "$q_file" "$x_file" "$out0" 0 "$LEAF_CAPACITY" "$REL_ERROR" "$B_VALUE" > "$log0"
    "$BIN" "$q_file" "$x_file" "$out3" 3 "$LEAF_CAPACITY" "$REL_ERROR" "$B_VALUE" > "$log3"

    local m0_qps
    local m3_qps
    local m3_over_m0
    local mean_q
    local mean_x

    m0_qps="$(awk '/^Method 0:/{print $3}' "$log0" | tail -n 1)"
    m3_qps="$(awk '/^Method 3:/{print $3}' "$log3" | tail -n 1)"
    m3_over_m0="$(awk -v a="$m3_qps" -v b="$m0_qps" 'BEGIN{if(b==0) print "nan"; else printf "%.6f", a/b}')"
    mean_q="$(mean_abs_file "$q_file")"
    mean_x="$(mean_abs_file "$x_file")"

    echo "${dataset},${scale},${m0_qps},${m3_qps},${m3_over_m0},${mean_q},${mean_x}" >> "$SUMMARY_CSV"
}

for s in $SCALES; do
    s_tag="${s//./p}"

    HT_Q_S="/tmp/ht_q300_s${s_tag}.data"
    HT_X_S="/tmp/ht_x30000_s${s_tag}.data"
    MINI_Q_S="/tmp/miniboone_q300_s${s_tag}.data"
    MINI_X_S="/tmp/miniboone_x30000_s${s_tag}.data"

    make_scaled_file "$HT_Q" "$HT_Q_S" "$s"
    make_scaled_file "$HT_X" "$HT_X_S" "$s"
    make_scaled_file "$MINI_Q" "$MINI_Q_S" "$s"
    make_scaled_file "$MINI_X" "$MINI_X_S" "$s"

    echo "[RUN] HT scale=${s}"
    run_one "HT" "$HT_Q_S" "$HT_X_S" "$s"

    echo "[RUN] MiniBooNE scale=${s}"
    run_one "MiniBooNE" "$MINI_Q_S" "$MINI_X_S" "$s"
done

echo
echo "Summary CSV: $SUMMARY_CSV"
cat "$SUMMARY_CSV"
