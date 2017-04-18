#!/bin/bash

#replace "SEQHEADER" with the first 4 characters of your sequence name including the @
#replace "DELIMITER" with the specific delimiter from your fastq file, 
#see the two examples below:
#@HWI-ST141_0363:2:1101:1175:2080#ATCACG/1
#@HWI-ST1117:83:D0TKCACXX:6:1101:1274:2051 2:N:0:TAGCTT

#SEQHEADER would be "@HWI" for both examples, DELIMITER would be a "/" in the first case and
#a " " (space) in the second. Both SEQHEADER and DELIMITER need to be in "", while the filenames
#do not.

#Then copy the following line as many times as the number of paired end lanes that you have, then
#replace YOURFILE#1_1_trimmed_clipped.fastq with the name of your forward direction .fastq file
#replace YOURFILE#1_2_trimmed_clipped.fastq with the name of your reverse direction .fastq file 
#for that same sample.

SEQHEADER="@HIS"
DELIMITER=" "

# ========================== NEG1 ====================================== #

FILENAME_IN="neg1"
NUM_FILES_PER_READ=4

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastqcombinepairedend.py "@HIS" " " $FILENAME_IN"_R1_00"$i"_clipped.fastq" $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

# ========================== POS1 ====================================== #

FILENAME_IN="pos1"
NUM_FILES_PER_READ=5

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastqcombinepairedend.py "@HIS" " " $FILENAME_IN"_R1_00"$i"_clipped.fastq" $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

# ========================== POS2 ====================================== #

FILENAME_IN="pos2"
NUM_FILES_PER_READ=2

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastqcombinepairedend.py "@HIS" " " $FILENAME_IN"_R1_00"$i"_clipped.fastq" $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

##### END OF CONTROLS ... SAMPLES NEXT #####

# ========================== BAD1 ====================================== #

FILENAME_IN="bad1"
NUM_FILES_PER_READ=4

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastqcombinepairedend.py "@HIS" " " $FILENAME_IN"_R1_00"$i"_clipped.fastq" $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

# ========================== BAD2 ====================================== #

FILENAME_IN="bad2"
NUM_FILES_PER_READ=6

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastqcombinepairedend.py "@HIS" " " $FILENAME_IN"_R1_00"$i"_clipped.fastq" $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

# ========================== GOOD1 ====================================== #

FILENAME_IN="good1"
NUM_FILES_PER_READ=6

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastqcombinepairedend.py "@HIS" " " $FILENAME_IN"_R1_00"$i"_clipped.fastq" $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

# ========================== GOOD2 ====================================== #

FILENAME_IN="good2"
NUM_FILES_PER_READ=6

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastqcombinepairedend.py "@HIS" " " $FILENAME_IN"_R1_00"$i"_clipped.fastq" $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

# ========================== MIXED1 ====================================== #

FILENAME_IN="mixed1"
NUM_FILES_PER_READ=6

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastqcombinepairedend.py "@HIS" " " $FILENAME_IN"_R1_00"$i"_clipped.fastq" $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

# ========================== MIXED2 ====================================== #

FILENAME_IN="mixed2"
NUM_FILES_PER_READ=6

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastqcombinepairedend.py "@HIS" " " $FILENAME_IN"_R1_00"$i"_clipped.fastq" $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

