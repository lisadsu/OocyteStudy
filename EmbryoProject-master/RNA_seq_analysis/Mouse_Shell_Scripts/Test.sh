#!/bin/bash

# batch processing of all reads for a given folder
# define variables
FILE_TO_ANALYZE="JH_neg_IL4321-5_ACAGTG_L001"
FILENAME_OUT="neg1"
NUM_FILES_PER_READ=4

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILE_TO_ANALYZE"_R1_00"$i".fastq"
	echo $FILE_TO_ANALYZE"_R2_00"$i".fastq"
	echo $FILENAME_OUT"_R1_00"$i"_trimmed.fastq"
	echo $FILENAME_OUT"_R2_00"$i"_trimmed.fastq"

#	fastq_quality_trimmer 
	echo -v -t 20 -l 20 -i $FILE_TO_ANALYZE"_R1_00"$i".fastq" -o $FILENAME_OUT"_R1_00"$i"_trimmed.fastq"

#	fastq_quality_trimmer 
	echo -v -t 20 -l 20 -i $FILE_TO_ANALYZE"_R2_00"$i".fastq" -o $FILENAME_OUT"_R2_00"$i"_trimmed.fastq"

done

# echo ${fileNameIn_R1[@]}
# echo ${fileNameIn_R2[@]}
# echo ${fileNameOut_R1[@]}
# echo ${fileNameOut_R2[@]}







