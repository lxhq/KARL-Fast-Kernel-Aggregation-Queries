#!/usr/bin/env bash
set -euo pipefail

# Config
LEAF_CAPACITY="${LEAF_CAPACITY:-40}"
B_VALUES="${B_VALUES:-1 5 10}"
REL_ERRORS="${REL_ERRORS:-0.1 0.2 0.5}"
OUTPUT_ROOT="${OUTPUT_ROOT:-./test_Result/pruning_conditions}"

BIN="./build/release/main_release"
if [[ ! -x "$BIN" ]]; then
    echo "Error: release binary not found: $BIN" >&2
    exit 1
fi

mkdir -p "$OUTPUT_ROOT"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${OUTPUT_ROOT}/run_${RUN_ID}"
mkdir -p "$RUN_DIR"
SUMMARY_CSV="${RUN_DIR}/summary.csv"

# Create local subsets (required due runtime sandbox path restrictions).
make_subset() {
    local q_src="$1"
    local x_src="$2"
    local q_out="$3"
    local x_out="$4"
    local q_num="$5"
    local x_num="$6"
    local dim="$7"

    { echo "${q_num} ${dim}"; sed -n "2,$((q_num+1))p" "$q_src"; } > "$q_out"
    { echo "${x_num} ${dim}"; sed -n "2,$((x_num+1))p" "$x_src"; } > "$x_out"
}

Q_HT="/tmp/karl_exp_ht_q300.data"
X_HT="/tmp/karl_exp_ht_x30000.data"
Q_SUSY="/tmp/karl_exp_susy_q300.data"
X_SUSY="/tmp/karl_exp_susy_x30000.data"
Q_MINI="/tmp/karl_exp_miniboone_q300.data"
X_MINI="/tmp/karl_exp_miniboone_x30000.data"
Q_MNIST="/tmp/karl_exp_mnist_q100.data"
X_MNIST="/tmp/karl_exp_mnist_x10000.data"

make_subset \
    "/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/home/HT_Sensor_dataset_sampleSet.data" \
    "/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/home/HT_Sensor_dataset_X.data" \
    "$Q_HT" "$X_HT" 300 30000 10

make_subset \
    "/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/susy/SUSY_sampleSet.data" \
    "/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/susy/SUSY_X.data" \
    "$Q_SUSY" "$X_SUSY" 300 30000 18

make_subset \
    "/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/miniboone/MiniBooNE_sampleSet.data" \
    "/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/miniboone/MiniBooNE_X.data" \
    "$Q_MINI" "$X_MINI" 300 30000 50

make_subset \
    "/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/mnist/MNIST_qSet.data" \
    "/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/mnist/MNIST_X.data" \
    "$Q_MNIST" "$X_MNIST" 100 10000 784

echo "dataset,dim,q_num,n_num,b,rel_error,leaf,m0_qps,m3_qps,m3_over_m0" > "$SUMMARY_CSV"

run_case() {
    local dataset="$1"
    local dim="$2"
    local q_num="$3"
    local n_num="$4"
    local q_file="$5"
    local x_file="$6"
    local b="$7"
    local rel="$8"

    local key="${dataset}_b${b}_e${rel}"
    local out0="${RUN_DIR}/${key}_m0.txt"
    local out3="${RUN_DIR}/${key}_m3.txt"
    local log0="${RUN_DIR}/${key}_m0.log"
    local log3="${RUN_DIR}/${key}_m3.log"

    echo "[RUN] dataset=${dataset} dim=${dim} b=${b} rel_error=${rel}"
    "$BIN" "$q_file" "$x_file" "$out0" 0 "$LEAF_CAPACITY" "$rel" "$b" > "$log0"
    "$BIN" "$q_file" "$x_file" "$out3" 3 "$LEAF_CAPACITY" "$rel" "$b" > "$log3"

    local m0_qps
    local m3_qps
    local ratio

    m0_qps="$(awk '/^Method 0:/{print $3}' "$log0" | tail -n 1)"
    m3_qps="$(awk '/^Method 3:/{print $3}' "$log3" | tail -n 1)"
    ratio="$(awk -v a="$m3_qps" -v b="$m0_qps" 'BEGIN{if(b==0) print "nan"; else printf "%.6f", a/b}')"

    echo "${dataset},${dim},${q_num},${n_num},${b},${rel},${LEAF_CAPACITY},${m0_qps},${m3_qps},${ratio}" >> "$SUMMARY_CSV"
}

for rel in $REL_ERRORS; do
    for b in $B_VALUES; do
        run_case "HT" 10 300 30000 "$Q_HT" "$X_HT" "$b" "$rel"
        run_case "SUSY" 18 300 30000 "$Q_SUSY" "$X_SUSY" "$b" "$rel"
        run_case "MiniBooNE" 50 300 30000 "$Q_MINI" "$X_MINI" "$b" "$rel"
        run_case "MNIST" 784 100 10000 "$Q_MNIST" "$X_MNIST" "$b" "$rel"
    done
done

echo
echo "Summary CSV: $SUMMARY_CSV"
cat "$SUMMARY_CSV"
