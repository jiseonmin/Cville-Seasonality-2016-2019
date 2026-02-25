#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --output=logs/sra_download_%A_%a.log
#SBATCH --error=logs/sra_download_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=lotterhos

# ── User-defined paths ────────────────────────────────────────────────────────
SRR_LIST="../../DATA/individual_SRA.txt"
METADATA_CSV="../../DATA/Table_S1_GENETICS-2023-306673.csv"
REF_FILE="./Guide_Files/Alyssa_ind_reads_guideFile.txt"
SRA_DIR="../../DATA/sra_files"
FASTQ_DIR="../../DATA/fastq_raw"
FINAL_DIR="../../DATA/fastq_renamed"
THREADS=${SLURM_CPUS_PER_TASK}
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p logs "$SRA_DIR" "$FASTQ_DIR" "$FINAL_DIR"

module load sratoolkit

# ── Get this job's SRR accession ─────────────────────────────────────────────
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST" | tr -d '[:space:]')

echo "=========================================="
echo "Job started: $(date)"
echo "Array task:  $SLURM_ARRAY_TASK_ID"
echo "SRR:         $SRR"
echo "Threads:     $THREADS"
echo "=========================================="

# ── Step 1: Build SRR → sampleId lookup from metadata CSV ────────────────────
# Only loads rows where sampleId starts with "CM_"
declare -A SRR_TO_SAMPLE
while IFS=',' read -r sampleId country city locality collectionDate nFlies SRA_accession rest; do
    [[ "$sampleId" == "sampleId" ]] && continue   # skip header
    [[ -z "$SRA_accession" ]] && continue
    [[ "$sampleId" != CM_* ]] && continue          # only CM_ samples
    SRR_TO_SAMPLE["$SRA_accession"]="$sampleId"
done < "$METADATA_CSV"

# ── Step 2: Build sampleId (dot format) → HL prefix lookup from ref file ─────
declare -A DOT_SAMPLE_TO_R1
declare -A DOT_SAMPLE_TO_R2

while IFS=$'\t' read -r r1_full lane1 r1_base r2_full lane2 r2_base dot_sample sl_suffix no_ext; do
    [[ -z "$dot_sample" ]] && continue
    DOT_SAMPLE_TO_R1["$dot_sample"]="$r1_full"
    DOT_SAMPLE_TO_R2["$dot_sample"]="$r2_full"
done < "$REF_FILE"

# ── Helper: convert sampleId underscores → dots ───────────────────────────────
underscore_to_dot() {
    echo "${1//_/.}"
}

# ── Step 3: prefetch ──────────────────────────────────────────────────────────
echo ""
echo "--- Step 3: prefetch ---"
SRA_PATH="$SRA_DIR/$SRR/$SRR.sra"
if [[ -f "$SRA_PATH" ]]; then
    echo "  [SKIP] $SRR already downloaded."
else
    echo "  Fetching $SRR ..."
    prefetch "$SRR" --output-directory "$SRA_DIR" --max-size 100GB
    if [[ $? -ne 0 ]]; then
        echo "  [ERROR] prefetch failed for $SRR" >&2
        exit 1
    fi
fi

# ── Step 4: fasterq-dump + gzip + rename ─────────────────────────────────────
echo ""
echo "--- Step 4: fasterq-dump + gzip + rename ---"

SAMPLE="${SRR_TO_SAMPLE[$SRR]}"
if [[ -z "$SAMPLE" ]]; then
    echo "  [WARN] No sampleId found for $SRR — using SRR name as fallback."
    SAMPLE="$SRR"
fi

DOT_SAMPLE=$(underscore_to_dot "$SAMPLE")
R1_TARGET="${DOT_SAMPLE_TO_R1[$DOT_SAMPLE]}"
R2_TARGET="${DOT_SAMPLE_TO_R2[$DOT_SAMPLE]}"

# Check if final renamed files already exist
if [[ -n "$R1_TARGET" && -f "$FINAL_DIR/$R1_TARGET" ]]; then
    echo "  [SKIP] $SRR already processed → $R1_TARGET"
    exit 0
fi

echo "  Processing $SRR (sample: $SAMPLE, dot: $DOT_SAMPLE) ..."

fasterq-dump "$SRA_DIR/$SRR/$SRR.sra" \
    --outdir "$FASTQ_DIR" \
    --threads "$THREADS" \
    --split-files \
    --skip-technical

if [[ $? -ne 0 ]]; then
    echo "  [ERROR] fasterq-dump failed for $SRR" >&2
    exit 1
fi

for READ in 1 2; do
    RAW_FQ="$FASTQ_DIR/${SRR}_${READ}.fastq"
    if [[ ! -f "$RAW_FQ" ]]; then
        echo "  [WARN] Expected $RAW_FQ not found — may be single-end." >&2
        continue
    fi

    if [[ "$READ" -eq 1 && -n "$R1_TARGET" ]]; then
        TARGET_NAME="$R1_TARGET"
    elif [[ "$READ" -eq 2 && -n "$R2_TARGET" ]]; then
        TARGET_NAME="$R2_TARGET"
    else
        TARGET_NAME="${SRR}_${READ}.fastq.gz"
        echo "  [WARN] No HL reference name for $DOT_SAMPLE R${READ} — using $TARGET_NAME"
    fi

    echo "    Compressing → $FINAL_DIR/$TARGET_NAME"
    gzip -c "$RAW_FQ" > "$FINAL_DIR/$TARGET_NAME"
    rm -f "$RAW_FQ"
done

echo ""
echo "=========================================="
echo "Job finished: $(date)"
echo "=========================================="
