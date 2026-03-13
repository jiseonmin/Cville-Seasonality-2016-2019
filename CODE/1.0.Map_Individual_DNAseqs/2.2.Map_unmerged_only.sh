#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=./logs/Step2_Map_Unmerged_%A_%a.out
#SBATCH --error=./logs/Step2_Map_Unmerged_%A_%a.err
#SBATCH --partition=lotterhos

# STEP 2: Map unmerged reads (BWA -> samtools -> picard SortSam -> MarkDuplicates)

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

echo "Processing unmerged reads for sample: ${i}"

echo "Running BWA for unmerged reads..."
bwa mem -M -t $CPU \
	-R "@RG\tID:${i}\tSM:${i}\tLB:${i}\tPL:ILLUMINA" \
	$REFERENCE \
    $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.trim.1.fq \
    $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.trim.2.fq \
    > $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.sam

echo "Running samtools flagstat..."
samtools flagstat --threads $CPU \
    $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.sam \
    > $WORKING_FOLDER/unmerged_reads/${i}/${i}.flagstats_raw_unmerged.sam.txt

echo "Running samtools view..."
samtools view -b -q $QUAL --threads $CPU \
    $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.sam \
    > $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.bam

echo "Running Picard SortSam..."
java -Xmx$JAVAMEM -jar $PICARD SortSam \
    I=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.bam \
    O=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.srt.bam \
    SO=coordinate VALIDATION_STRINGENCY=SILENT

echo "Running Picard MarkDuplicates..."
java -Xmx$JAVAMEM -jar $PICARD MarkDuplicates \
    I=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.srt.bam \
    O=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.srt.rmdp.bam \
    M=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.dupstat.txt \
    VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

echo "Running Qualimap..."
qualimap bamqc \
    -bam $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.srt.rmdp.bam \
    -outdir $WORKING_FOLDER/mapping_stats/Qualimap_${i}_unmerged \
    --java-mem-size=$JAVAMEM

echo "Cleaning up intermediates..."
rm $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.sam
rm $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.bam
rm $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.srt.bam

mv $WORKING_FOLDER/unmerged_reads/${i}/${i}.flagstats_raw_unmerged.sam.txt \
    $WORKING_FOLDER/mapping_stats
mv $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.dupstat.txt \
    $WORKING_FOLDER/mapping_stats

echo "Step 2 complete for ${i}"
