#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=8:00:00
#SBATCH --output=./logs/Step1_Trim_%A_%a.out
#SBATCH --error=./logs/Step1_Trim_%A_%a.err
#SBATCH --partition=lotterhos

# STEP 1: Trim merged and unmerged reads with BBDuk

module load bbmap

WORKING_FOLDER=/scratch/j.min/data/fastq-to-vcf
SAMPLE_FILE="/projects/lotterhos/Cville-Seasonality-2016-2019/CODE/1.0.Map_Individual_DNAseqs/Guide_Files/Alyssa_ind_reads_guideFile.txt"

i=`awk -F "\t" '{print $9}' $SAMPLE_FILE | sed -n ${SLURM_ARRAY_TASK_ID}p`

echo "Trimming merged reads for ${i}"
bbduk.sh \
    in=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq \
    out=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.trim.fq \
    ftl=15 ftr=285 qtrim=w trimq=20

rm $WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq

echo "Trimming unmerged reads for ${i}"
bbduk.sh \
    in=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.1.fq \
    in2=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.2.fq \
    out=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.trim.1.fq \
    out2=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.trim.2.fq \
    ftl=15 qtrim=w trimq=20

rm $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.1.fq
rm $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.2.fq

echo "Step 1 complete for ${i}"
