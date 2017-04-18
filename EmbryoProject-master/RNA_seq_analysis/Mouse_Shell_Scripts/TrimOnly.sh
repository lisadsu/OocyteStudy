#!/bin/bash
#This script is a multi-script that quality trims and then clips adapter seqs from .fastq files
#Quality trimming:
#Quality score cutoff 20, then remove any sequences that are shorter than 20 bases after trimming.

# batch processing of all reads for a given folder
# Steps:
# 	1. Fill in values for variables below (do this for each new folder of data)
# 	2. Navigate to appropriate directory where data is stored


# define variables
FILE_TO_ANALYZE="JH_neg_IL4321-5_ACAGTG_L001"
FILENAME_OUT="neg1"
NUM_FILES_PER_READ=4

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILE_TO_ANALYZE"_R1_00"$i".fastq"
	echo $FILENAME_OUT"_R1_00"$i"_trimmed.fastq"

	fastq_quality_trimmer -v -t 20 -l 20 -i $FILE_TO_ANALYZE"_R1_00"$i".fastq" -o $FILENAME_OUT"_R1_00"$i"_trimmed.fastq"

	echo $FILE_TO_ANALYZE"_R2_00"$i".fastq"
	echo $FILENAME_OUT"_R2_00"$i"_trimmed.fastq"

	fastq_quality_trimmer -v -t 20 -l 20 -i $FILE_TO_ANALYZE"_R2_00"$i".fastq" -o $FILENAME_OUT"_R2_00"$i"_trimmed.fastq"
done

mkdir ./Trimmed
mv *trimmed.fastq ./Trimmed
cd ./Trimmed

