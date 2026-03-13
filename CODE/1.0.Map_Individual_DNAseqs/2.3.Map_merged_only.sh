#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=./logs/Step3_Map_Merged_%A_%a.out
#SBATCH --error=./logs/Step3_Map_Merged_%A_%a.err
#SBATCH --partition=lotterhos

# STEP 3: Map merged reads + MergeSamFiles

set -e

module load bwa
module load samtools
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate cville

PICARD=$(find $CONDA_PREFIX -name "picard.jar" | head -1)

WORKING_FOLDER=/scratch/j.min/data/fastq-to-vcf
SAMPLE_FILE="/projects/lotterhos/Cville-Seasonality-2016-2019/CODE/1.0.Map_Individual_DNAseqs/Guide_Files/Alyssa_ind_reads_guideFile.txt"
REFERENCE=/scratch/j.min/data/holo_dmel_6.12.fa
JAVAMEM=35g
CPU=$SLURM_CPUS_ON_NODE
QUAL=40

i=`awk -F "\t" '{print $9}' $SAMPLE_FILE | sed -n ${SLURM_ARRAY_TASK_ID}p`

echo "Processing sample: ${i}"

# ── BWA + samtools + picard for merged (skip if already done) ────────────────
if [[ -f "$WORKING_FOLDER/merged_reads/${i}/${i}.merged.srt.rmdp.bam" ]]; then
    echo "merged.srt.rmdp.bam already exists — skipping merged mapping steps."
else
    echo "Running BWA for merged reads..."
    bwa mem -M -t $CPU \
        -R "@RG\tID:${i}\tSM:${i}\tLB:${i}\tPL:ILLUMINA" \
        $REFERENCE \
        $WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.trim.fq \
        > $WORKING_FOLDER/merged_reads/${i}/${i}.merged.sam

    samtools flagstat --threads $CPU \
        $WORKING_FOLDER/merged_reads/${i}/${i}.merged.sam \
        > $WORKING_FOLDER/merged_reads/${i}/${i}.flagstats_raw_merged.sam.txt

    samtools view -b -q $QUAL --threads $CPU \
        $WORKING_FOLDER/merged_reads/${i}/${i}.merged.sam \
        > $WORKING_FOLDER/merged_reads/${i}/${i}.merged.bam

    java -Xmx$JAVAMEM -jar $PICARD SortSam \
        I=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.bam \
        O=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.srt.bam \
        SO=coordinate VALIDATION_STRINGENCY=SILENT

    java -Xmx$JAVAMEM -jar $PICARD MarkDuplicates \
        I=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.srt.bam \
        O=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.srt.rmdp.bam \
        M=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.dupstat.txt \
        VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

    rm $WORKING_FOLDER/merged_reads/${i}/${i}.merged.sam
    rm $WORKING_FOLDER/merged_reads/${i}/${i}.merged.bam
    rm $WORKING_FOLDER/merged_reads/${i}/${i}.merged.srt.bam

    mv $WORKING_FOLDER/merged_reads/${i}/${i}.flagstats_raw_merged.sam.txt \
        $WORKING_FOLDER/mapping_stats
    mv $WORKING_FOLDER/merged_reads/${i}/${i}.merged.dupstat.txt \
        $WORKING_FOLDER/mapping_stats
fi

# ── MergeSamFiles (requires unmerged.srt.rmdp.bam from step 3) ───────────────
echo "Checking that unmerged.srt.rmdp.bam exists before merging..."
if [[ ! -f "$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.srt.rmdp.bam" ]]; then
    echo "ERROR: unmerged.srt.rmdp.bam not found. Run step2 first!" >&2
    exit 1
fi

echo "Running MergeSamFiles..."
java -Xmx$JAVAMEM -jar $PICARD MergeSamFiles \
    I=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.srt.rmdp.bam \
    I=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.srt.rmdp.bam \
    O=$WORKING_FOLDER/joint_bams/${i}.joint.bam

echo "Running SortSam on joint bam..."
java -Xmx$JAVAMEM -jar $PICARD SortSam \
    I=$WORKING_FOLDER/joint_bams/${i}.joint.bam \
    O=$WORKING_FOLDER/joint_bams/${i}.joint.srt.bam \
    SO=coordinate VALIDATION_STRINGENCY=SILENT

echo "Running MarkDuplicates on joint bam..."
java -Xmx$JAVAMEM -jar $PICARD MarkDuplicates \
    I=$WORKING_FOLDER/joint_bams/${i}.joint.srt.bam \
    O=$WORKING_FOLDER/joint_bams/${i}.joint.srt.rmdp.bam \
    M=$WORKING_FOLDER/mapping_stats/${i}.joint.dupstat.txt \
    VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

echo "Running Qualimap on joint bam..."
qualimap bamqc \
    -bam $WORKING_FOLDER/joint_bams/${i}.joint.srt.rmdp.bam \
    -outdir $WORKING_FOLDER/joint_bams_qualimap/Qualimap_JointBam_${i} \
    --java-mem-size=$JAVAMEM

rm $WORKING_FOLDER/joint_bams/${i}.joint.bam
rm $WORKING_FOLDER/joint_bams/${i}.joint.srt.bam

echo "Step 3 complete for ${i}"
