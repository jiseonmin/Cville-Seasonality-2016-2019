#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --output=logs/sra_download_%j.log
#SBATCH --error=logs/sra_download_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=lotterhos

# ── User-defined paths ────────────────────────────────────────────────────────
SRR_LIST="../../DATA/individual_SRA.txt"                            # one SRR accession per line
METADATA_TSV="../../DATA/Table_S1_GENETICS-2023-306673.csv"         # the full table (sampleId + SRA_accession cols)
REF_FILE="./Guide_Files/Alyssa_ind_reads_guideFile.txt"            # the HL... reference file from the pipeline
SRA_DIR="../../DATA/sra_files"                  # where prefetch stores .sra files
FASTQ_DIR="../../DATA/fastq_raw"                # fasterq-dump output
FINAL_DIR="../../DATA/fastq_renamed"            # gzipped + renamed output
THREADS=${SLURM_CPUS_PER_TASK}
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p logs "$SRA_DIR" "$FASTQ_DIR" "$FINAL_DIR"

# Load modules — adjust module names to match your HPC
module load sratoolkit

echo "=========================================="
echo "Job started: $(date)"
echo "SRR list: $SRR_LIST"
echo "Threads: $THREADS"
echo "=========================================="

# ── Step 1: Build SRR → sampleId lookup from metadata CSV ───────────────────
# Expects columns: sampleId (col 1) and SRA_accession (col 7) — adjust if needed
# Only loads rows where sampleId starts with "CM_"
declare -A SRR_TO_SAMPLE
while IFS=',' read -r sampleId country city locality collectionDate nFlies SRA_accession rest; do
    [[ "$sampleId" == "sampleId" ]] && continue   # skip header
    [[ -z "$SRA_accession" ]] && continue
    [[ "$sampleId" != CM_* ]] && continue          # only CM_ samples
    SRR_TO_SAMPLE["$SRA_accession"]="$sampleId"
done < "$METADATA_TSV"

echo "Loaded ${#SRR_TO_SAMPLE[@]} SRR→sampleId mappings from metadata."

# ── Step 2: Build sampleId (dot format) → HL prefix lookup from ref file ─────
# Reference file columns (tab-separated):
#   1: HL..._1_CM.xxx.fastq.gz   (R1 full name)
#   2: HL..._s#                  (flowcell+lane)
#   3: CM.xxx_SLyyy.fastq.gz     (R1 basename, dot-format sample)
#   4: HL..._2_CM.xxx.fastq.gz   (R2 full name)
#   5: HL..._s#
#   6: CM.xxx_SLyyy.fastq.gz     (R2 basename)
#   7: CM.xxx                    (dot-format sampleId only)
#   8: SLyyy.fastq.gz
#   9: HL..._CM.xxx_SLyyy        (no read number, no extension)
declare -A DOT_SAMPLE_TO_R1   # dot_sampleId → full R1 filename
declare -A DOT_SAMPLE_TO_R2   # dot_sampleId → full R2 filename

while IFS=$'\t' read -r r1_full lane1 r1_base r2_full lane2 r2_base dot_sample sl_suffix no_ext; do
    [[ -z "$dot_sample" ]] && continue
    DOT_SAMPLE_TO_R1["$dot_sample"]="$r1_full"
    DOT_SAMPLE_TO_R2["$dot_sample"]="$r2_full"
done < "$REF_FILE"

echo "Loaded ${#DOT_SAMPLE_TO_R1[@]} sample→filename mappings from reference file."

# ── Helper: convert sampleId underscores → dots (e.g. CM_002_0916 → CM.002.0916)
underscore_to_dot() {
    echo "${1//_/.}"
}

# ── Step 3: prefetch all SRRs ─────────────────────────────────────────────────
echo ""
echo "--- Step 3: prefetch ---"
while read -r SRR; do
    [[ -z "$SRR" || "$SRR" =~ ^# ]] && continue

    SRA_PATH="$SRA_DIR/$SRR/$SRR.sra"
    if [[ -f "$SRA_PATH" ]]; then
        echo "  [SKIP] $SRR already downloaded."
        continue
    fi

    echo "  Fetching $SRR ..."
    prefetch "$SRR" --output-directory "$SRA_DIR" --max-size 100GB
    if [[ $? -ne 0 ]]; then
        echo "  [ERROR] prefetch failed for $SRR" >&2
    fi
done < "$SRR_LIST"

# ── Step 4: fasterq-dump → gzip → rename ─────────────────────────────────────
echo ""
echo "--- Step 4: fasterq-dump + gzip + rename ---"
while read -r SRR; do
    [[ -z "$SRR" || "$SRR" =~ ^# ]] && continue

    SAMPLE="${SRR_TO_SAMPLE[$SRR]}"
    if [[ -z "$SAMPLE" ]]; then
        echo "  [WARN] No sampleId found for $SRR — skipping rename, keeping SRR name."
        SAMPLE="$SRR"
    fi

    DOT_SAMPLE=$(underscore_to_dot "$SAMPLE")
    R1_TARGET="${DOT_SAMPLE_TO_R1[$DOT_SAMPLE]}"
    R2_TARGET="${DOT_SAMPLE_TO_R2[$DOT_SAMPLE]}"

    # Check if final renamed files already exist
    if [[ -n "$R1_TARGET" && -f "$FINAL_DIR/${R1_TARGET}.gz" ]]; then
        echo "  [SKIP] $SRR already processed → $R1_TARGET.gz"
        continue
    fi

    echo "  Processing $SRR (sample: $SAMPLE, dot: $DOT_SAMPLE) ..."

    # fasterq-dump
    fasterq-dump "$SRA_DIR/$SRR/$SRR.sra" \
        --outdir "$FASTQ_DIR" \
        --threads "$THREADS" \
        --split-files \
        --skip-technical

    if [[ $? -ne 0 ]]; then
        echo "  [ERROR] fasterq-dump failed for $SRR" >&2
        continue
    fi

    # Compress with pigz (parallel), fall back to gzip
    for READ in 1 2; do
        RAW_FQ="$FASTQ_DIR/${SRR}_${READ}.fastq"
        if [[ ! -f "$RAW_FQ" ]]; then
            echo "  [WARN] Expected $RAW_FQ not found — may be single-end." >&2
            continue
        fi

        # Determine target filename
        if [[ "$READ" -eq 1 && -n "$R1_TARGET" ]]; then
            TARGET_NAME="${R1_TARGET}.gz"
        elif [[ "$READ" -eq 2 && -n "$R2_TARGET" ]]; then
            TARGET_NAME="${R2_TARGET}.gz"
        else
            # Fallback: no ref entry found, keep SRR-based name
            TARGET_NAME="${SRR}_${READ}.fastq.gz"
            echo "  [WARN] No HL reference name for $DOT_SAMPLE R${READ} — using $TARGET_NAME"
        fi

        echo "    Compressing → $FINAL_DIR/$TARGET_NAME"
        gzip -c "$RAW_FQ" > "$FINAL_DIR/$TARGET_NAME"

        rm -f "$RAW_FQ"   # remove uncompressed copy to save space
    done

done < "$SRR_LIST"

echo ""
echo "=========================================="
echo "All done: $(date)"
echo "Final files in: $FINAL_DIR"
echo "=========================================="
