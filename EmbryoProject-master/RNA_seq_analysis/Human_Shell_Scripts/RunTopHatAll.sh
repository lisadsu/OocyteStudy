#!/bin/bash

# Part 2 in RNA-seq analysis

# Separate out true paired reads from singled reads (fastqpairedendcombiner.py)
# Run TopHat on all files

BASE_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanMiSeq"
RAW_DATA_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanMiSeq/ClippedTrimmed"

declare -a sampleNum=(1 2 3 4 5 6 7 8 9 10 11 12)
declare -a embryoNum=(1 2 3 4 5 6 7 9 10 12 13 14)

#declare -a sampleNum=(1)
#declare -a embryoNum=(7)

RVAL=150

# for paired end combiner
SEQHEADER="@MIS"
DELIMITER=" "

cd $BASE_DIR

# uncomment this if no transcriptome index previously built by TopHat
# tophat2 -G $REF_GTF_HUMAN --transcriptome-index ./transcriptome_data/known $HUMAN_GENOME

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

		echo "not combining paired ends"

#		fastqcombinepairedend.py $SEQHEADER " " "./ClippedTrimmed/"$FILENAME_IN"_R1_00"$j"_cliptrim.fastq" "./ClippedTrimmed/"$FILENAME_IN"_R2_00"$j"_cliptrim.fastq"


	done

	##############################
	# Run Tophat for each sample #
	##############################

	echo ""
	echo "now aligning embryo : E"$CURR_EMBRYO_NUM

	tophat2 -r $RVAL -G $REF_GTF_HUMAN -o "./TophatOut_"$FILENAME_IN"_paired" --transcriptome-index ./transcriptome_data/known --transcriptome-only $HUMAN_GENOME "./ClippedTrimmed/"$FILENAME_IN"_R1_stillpaired.fastq" "./ClippedTrimmed/"$FILENAME_IN"_R2_stillpaired.fastq"

	tophat2 -r $RVAL -G $REF_GTF_HUMAN -o "./TophatOut_"$FILENAME_IN"_singles" --transcriptome-index ./transcriptome_data/known --transcriptome-only $HUMAN_GENOME "./ClippedTrimmed/"$FILENAME_IN"_singles.fastq"


done






















