#!/bin/bash
#This script first collapses the trimmed and clipped fastq files into fasta files with only unique reads, the number of identical reads for each unique in parentheses:

#Copy the following line as many times as the number of samples you are working with, then
#replace "YOURFILE" with the names of your samples e.g. "1_Hot":

# ========================== NEG1 ====================================== #

FILENAME_IN="neg1"
NUM_FILES_PER_READ=4

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R1_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R1_00"$i"_clipped.fastq" -o $FILENAME_IN"_R1_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R1_00"$i"_collapsed.txt" $FILENAME_IN"_R1_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R2_00"$i"_clipped.fastq" -o $FILENAME_IN"_R2_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R2_00"$i"_collapsed.txt" $FILENAME_IN"_R2_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

done

# ========================== POS1 ====================================== #

FILENAME_IN="pos1"
NUM_FILES_PER_READ=5

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R1_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R1_00"$i"_clipped.fastq" -o $FILENAME_IN"_R1_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R1_00"$i"_collapsed.txt" $FILENAME_IN"_R1_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R2_00"$i"_clipped.fastq" -o $FILENAME_IN"_R2_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R2_00"$i"_collapsed.txt" $FILENAME_IN"_R2_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

done

# ========================== POS2 ====================================== #

FILENAME_IN="pos2"
NUM_FILES_PER_READ=2

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R1_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R1_00"$i"_clipped.fastq" -o $FILENAME_IN"_R1_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R1_00"$i"_collapsed.txt" $FILENAME_IN"_R1_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R2_00"$i"_clipped.fastq" -o $FILENAME_IN"_R2_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R2_00"$i"_collapsed.txt" $FILENAME_IN"_R2_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

done

##### END OF CONTROLS ... SAMPLES NEXT #####

# ========================== BAD1 ====================================== #

FILENAME_IN="bad1"
NUM_FILES_PER_READ=4

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R1_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R1_00"$i"_clipped.fastq" -o $FILENAME_IN"_R1_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R1_00"$i"_collapsed.txt" $FILENAME_IN"_R1_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R2_00"$i"_clipped.fastq" -o $FILENAME_IN"_R2_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R2_00"$i"_collapsed.txt" $FILENAME_IN"_R2_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

done

# ========================== BAD2 ====================================== #

FILENAME_IN="bad2"
NUM_FILES_PER_READ=6

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R1_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R1_00"$i"_clipped.fastq" -o $FILENAME_IN"_R1_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R1_00"$i"_collapsed.txt" $FILENAME_IN"_R1_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R2_00"$i"_clipped.fastq" -o $FILENAME_IN"_R2_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R2_00"$i"_collapsed.txt" $FILENAME_IN"_R2_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

done

# ========================== GOOD1 ====================================== #

FILENAME_IN="good1"
NUM_FILES_PER_READ=6

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R1_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R1_00"$i"_clipped.fastq" -o $FILENAME_IN"_R1_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R1_00"$i"_collapsed.txt" $FILENAME_IN"_R1_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R2_00"$i"_clipped.fastq" -o $FILENAME_IN"_R2_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R2_00"$i"_collapsed.txt" $FILENAME_IN"_R2_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

done

# ========================== GOOD2 ====================================== #

FILENAME_IN="good2"
NUM_FILES_PER_READ=6

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R1_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R1_00"$i"_clipped.fastq" -o $FILENAME_IN"_R1_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R1_00"$i"_collapsed.txt" $FILENAME_IN"_R1_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R2_00"$i"_clipped.fastq" -o $FILENAME_IN"_R2_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R2_00"$i"_collapsed.txt" $FILENAME_IN"_R2_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

done

# ========================== MIXED1 ====================================== #

FILENAME_IN="mixed1"
NUM_FILES_PER_READ=6

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R1_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R1_00"$i"_clipped.fastq" -o $FILENAME_IN"_R1_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R1_00"$i"_collapsed.txt" $FILENAME_IN"_R1_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R2_00"$i"_clipped.fastq" -o $FILENAME_IN"_R2_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R2_00"$i"_collapsed.txt" $FILENAME_IN"_R2_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

done

# ========================== MIXED2 ====================================== #

FILENAME_IN="mixed2"
NUM_FILES_PER_READ=6

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R1_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R1_00"$i"_clipped.fastq" -o $FILENAME_IN"_R1_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R1_00"$i"_collapsed.txt" $FILENAME_IN"_R1_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_collapsed.txt"

	fastx_collapser -v -i $FILENAME_IN"_R2_00"$i"_clipped.fastq" -o $FILENAME_IN"_R2_00"$i"_collapsed.txt"
	fastqduplicatecounter.py $FILENAME_IN"_R2_00"$i"_collapsed.txt" $FILENAME_IN"_R2_00"$i"_collapsed_headers.txt" >> $FILENAME_IN"_R1_00"$i"_duplicateCount.txt"

done


#This last line removes the dummy headerfiles, which are no longer necessary.
rm *_collapsed_headers.txt


# mkdir ../DuplicateCount
mv *duplicateCount.txt ../DuplicateCount


