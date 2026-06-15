#!/usr/bin/env bash
set -euo pipefail

# Set your test inputs here.
QUERY_FILE="/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/home/HT_Sensor_dataset_sampleSet.data"
DATASET_FILE="/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/home/HT_Sensor_dataset_X.data"
# QUERY_FILE="/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/miniboone/MiniBooNE_sampleSet.data"
# DATASET_FILE="/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/miniboone/MiniBooNE_X.data"
# QUERY_FILE="/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/susy/SUSY_sampleSet_100.data"
# DATASET_FILE="/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/susy/SUSY_X.data"
# QUERY_FILE="/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/mnist/MNIST_sampleSet.data"
# DATASET_FILE="/home/ubuntu/Documents/workspace/dataset/GPU-accelerated_Kernel_Density_Estimation/mnist/MNIST_X.data"
REL_ERROR="0.2"
OUTPUT_DIR="./test_Result/leaf_capacity"
METHODS=(0)
B_VALUES=(0.1 0.5 1 5 10)
LEAF_CAPACITIES=(10)
DATASET_BASENAME="$(basename "$DATASET_FILE")"
DATASET_NAME="${DATASET_BASENAME%.*}"

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
METHOD_FIRST="${METHODS[0]}"
METHOD_LAST="${METHODS[$((${#METHODS[@]}-1))]}"
B_FIRST="${B_VALUES[0]}"
B_LAST="${B_VALUES[$((${#B_VALUES[@]}-1))]}"
LEAF_FIRST="${LEAF_CAPACITIES[0]}"
LEAF_LAST="${LEAF_CAPACITIES[$((${#LEAF_CAPACITIES[@]}-1))]}"
LOG_FILE="${OUTPUT_DIR}/${DATASET_NAME}_method_${METHOD_FIRST}_to_${METHOD_LAST}_leaf_${LEAF_FIRST}_to_${LEAF_LAST}_eps_${REL_ERROR}_b_${B_FIRST}_to_${B_LAST}_${TIMESTAMP}.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "Binary: $BIN"
echo "Query file: $QUERY_FILE"
echo "Dataset file: $DATASET_FILE"
echo "Relative error: $REL_ERROR"
echo "Methods: ${METHODS[*]}"
echo "B values: ${B_VALUES[*]}"
echo "Output dir: $OUTPUT_DIR"
echo "Log file: $LOG_FILE"
echo

for METHOD in "${METHODS[@]}"; do
    for B_VALUE in "${B_VALUES[@]}"; do
        for LEAF in "${LEAF_CAPACITIES[@]}"; do
            OUT_KD="${OUTPUT_DIR}/${DATASET_NAME}_method_${METHOD}_leaf_${LEAF}_eps_${REL_ERROR}_b_${B_VALUE}.txt"

            echo "[KARL]   method=${METHOD}, leafCapacity=${LEAF}, b=${B_VALUE}"
            "$BIN" "$QUERY_FILE" "$DATASET_FILE" "$OUT_KD" "$METHOD" "$LEAF" "$REL_ERROR" "$B_VALUE"
        done
        echo
    done
done

echo
echo "Done. Results written to: $OUTPUT_DIR"
