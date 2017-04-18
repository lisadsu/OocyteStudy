#!/bin/bash

BASE_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/MouseHiSeq"
cd $BASE_DIR
#mkdir -p "./HTSeq_Count_Out"

declare -a sampleNum=("bad1" "bad2" "good1" "good2")

for SAMPLE in ${sampleNum[@]}
do

	FILENAME=$SAMPLE
	#mkdir -p "./"$FILENAME"_bam"

	echo "now sorting .bam file and running HTSeq-count for sample "$FILENAME

	# sort bam file
#	samtools sort -n "./Aligned/Bam/"$FILENAME".bam" "./"$FILENAME"_bam/"$FILENAME"_sorted"
	
	# then run HTSeq-count
	python -m HTSeq.scripts.count -f bam -o "./"$FILENAME"_bam/report_yes.sam" --stranded=yes "./"$FILENAME"_bam/"$FILENAME"_pairedonly_sorted.bam" $REF_GTF_MOUSE > "./HTSeq_Count_Out/"$FILENAME"_count_yes.txt"

	python -m HTSeq.scripts.count -f bam -o "./"$FILENAME"_bam/report_reverse.sam" --stranded=reverse "./"$FILENAME"_bam/"$FILENAME"_pairedonly_sorted.bam" $REF_GTF_MOUSE > "./HTSeq_Count_Out/"$FILENAME"_count_reverse.txt"

	python -m HTSeq.scripts.count -f bam -o "./"$FILENAME"_bam/report_no.sam" --stranded=no "./"$FILENAME"_bam/"$FILENAME"_pairedonly_sorted.bam" $REF_GTF_MOUSE > "./HTSeq_Count_Out/"$FILENAME"_count_no.txt"


done
