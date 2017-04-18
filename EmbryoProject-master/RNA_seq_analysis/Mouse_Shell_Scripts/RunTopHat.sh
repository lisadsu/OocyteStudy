#!/bin/bash

# Part 2 in RNA-seq analysis
# Run TopHat on all files

BASE_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/MouseHiSeq"
RAW_DATA_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/MouseHiSeq/ClippedTrimmed"

declare -a sampleNum=(1 2 3 4 5 6)
declare -a sampleName=("bad1" "bad2" "good1" "good2" "mixed1" "mixed2")

RVAL=150

cd $BASE_DIR

# uncomment this if no transcriptome index previously built by TopHat
# tophat2 -G $REF_GTF_HUMAN --transcriptome-index ./transcriptome_data/known $HUMAN_GENOME

# loop over samples
for i in ${sampleNum[@]}
do
	##############################
	# Run Tophat for each sample #
	##############################

	FILENAME_IN=${sampleName[$i - 1]}

	echo ""
	echo "now running TopHat for sample : "$FILENAME_IN

	tophat2 --no-discordant --no-mixed -p 4 -r $RVAL -G $REF_GTF_MOUSE -o "./TophatOut_"$FILENAME_IN"_paired" --transcriptome-index "./transcriptome_data/known" --transcriptome-only $MOUSE_GENOME "./ClippedTrimmed/"$FILENAME_IN"_R1_stillpaired.fastq" "./ClippedTrimmed/"$FILENAME_IN"_R2_stillpaired.fastq"

	tophat2 --no-discordant --no-mixed -p 4 -r $RVAL -G $REF_GTF_MOUSE -o "./TophatOut_"$FILENAME_IN"_singles" --transcriptome-index "./transcriptome_data/known" --transcriptome-only $MOUSE_GENOME "./ClippedTrimmed/"$FILENAME_IN"_singles.fastq"


done

