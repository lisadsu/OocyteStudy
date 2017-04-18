#!/bin/bash

# Part 2 in RNA-seq analysis

# Separate .fastq files into paired end reads and singles for better TopHat results

BASE_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanMiSeq"

declare -a sampleNum=(1 2 3 4 5 6 7 8 9 10 11 12)
declare -a embryoNum=(1 2 3 4 5 6 7 9 10 12 13 14)

# for paired end combiner
SEQHEADER="@MIS"
DELIMITER=" "

cd $BASE_DIR

# loop over samples
for i in ${sampleNum[@]}
do

	#######################################
	# Paired end combiner for each sample #
	#######################################

	CURR_EMBRYO_NUM=${embryoNum[$i - 1]}
	echo ""
	echo "now consolidating pairs for embryo : E"$CURR_EMBRYO_NUM

	NUM_FILES_PER_READ=1
	FILENAME_IN="E"$CURR_EMBRYO_NUM

	for ((j=1; j<=$NUM_FILES_PER_READ; j++))
	do

		fastqcombinepairedend.py "@MIS" " " "./ClippedTrimmed/"$FILENAME_IN"_R1_00"$j"_cliptrim.fastq" "./ClippedTrimmed/"$FILENAME_IN"_R2_00"$j"_cliptrim.fastq"

	done

done















