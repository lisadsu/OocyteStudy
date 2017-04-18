#!/bin/bash

#Adapter clipping:
#Matches to potential adapter sequences at the end of reads, then removes reads shorter than 20 bases. 

#For the next step, copy the line as many times as the number of samples, then 
#add the sequences of your adapters after the -a flag (e.g. -a ATTGGCTTTGGGCAT), as well as changing the sample names:

#mkdir ../Clipped

# ========================== NEG1 ====================================== #

FILENAME_IN="neg1"
NUM_FILES_PER_READ=4
ADAPTER_SEQUENCE="GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R1_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	echo $FILENAME_IN"_R2_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R2_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

# ========================== POS1 ====================================== #

FILENAME_IN="pos1"
NUM_FILES_PER_READ=5
ADAPTER_SEQUENCE="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCGATCTCGTATGCCGTCTTCTGCTTG"

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R1_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	echo $FILENAME_IN"_R2_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R2_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

# ========================== POS2 ====================================== #

FILENAME_IN="pos2"
NUM_FILES_PER_READ=2
ADAPTER_SEQUENCE="GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG"

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R1_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	echo $FILENAME_IN"_R2_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R2_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

##### END OF CONTROLS ... SAMPLES NEXT #####

# ========================== BAD1 ====================================== #

FILENAME_IN="bad1"
NUM_FILES_PER_READ=4
ADAPTER_SEQUENCE="GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTATGCCGTCTTCTGCTTG"

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R1_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	echo $FILENAME_IN"_R2_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R2_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

# ========================== BAD2 ====================================== #

FILENAME_IN="bad2"
NUM_FILES_PER_READ=6
ADAPTER_SEQUENCE="GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG"

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R1_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	echo $FILENAME_IN"_R2_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R2_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

# ========================== GOOD1 ====================================== #

FILENAME_IN="good1"
NUM_FILES_PER_READ=6
ADAPTER_SEQUENCE="GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R1_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	echo $FILENAME_IN"_R2_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R2_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

# ========================== GOOD2 ====================================== #

FILENAME_IN="good2"
NUM_FILES_PER_READ=6
ADAPTER_SEQUENCE="GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG"

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R1_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	echo $FILENAME_IN"_R2_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R2_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

# ========================== MIXED1 ====================================== #

FILENAME_IN="mixed1"
NUM_FILES_PER_READ=6
ADAPTER_SEQUENCE="GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCACATCTCGTATGCCGTCTTCTGCTTG"

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R1_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	echo $FILENAME_IN"_R2_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R2_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

# ========================== MIXED2 ====================================== #

FILENAME_IN="mixed2"
NUM_FILES_PER_READ=6
ADAPTER_SEQUENCE="GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG"

for ((i=1; i<=$NUM_FILES_PER_READ; i++))
do
	echo $FILENAME_IN"_R1_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R1_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R1_00"$i"_clipped.fastq"

	echo $FILENAME_IN"_R2_00"$i"_trimmed.fastq"
	echo $FILENAME_IN"_R2_00"$i"_clipped.fastq"

	fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i $FILENAME_IN"_R2_00"$i"_trimmed.fastq" -o $FILENAME_IN"_R2_00"$i"_clipped.fastq"
done

mv *clipped.fastq ../Clipped
cd ../Clipped




