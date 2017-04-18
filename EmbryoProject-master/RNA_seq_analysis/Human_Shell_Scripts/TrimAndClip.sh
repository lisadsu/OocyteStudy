#!/bin/bash

# Part 1 in RNA-seq analysis
# Trim ends with low quality and clip adapter sequences

RAW_DATA_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq"
PROCESSED_DATA_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/ClippedTrimmed"

declare -a sampleNum=(1 2 3 4 5 6 7 8 9 10 11 12)
declare -a embryoNum=(1 2 3 4 5 6 7 9 10 12 13 14)
declare -a shortAdaptArr=("ATCACG" "CGATGT" "TTAGGC" "TGACCA" "ACAGTG" "GCCAAT" "CAGATC" "ACTTGA" "GATCAG" "TAGCTT" "GGCTAC" "CTTGTA")
declare -a adaptArr=("GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG" "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG" "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG" "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG" "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG" "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG" "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG" "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG" "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG" "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG" "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG" "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG")

cd $RAW_DATA_DIR
mkdir -p ./ClippedTrimmed

# loop over samples
for i in ${sampleNum[@]}
do

	################
	# Part 1: trim #
	################

	CURR_EMBRYO_NUM=${embryoNum[$i - 1]}
	ADAPTER_SEQUENCE=${adaptArr[$i - 1]}

	echo ""
	echo "now trimming sample : E"$CURR_EMBRYO_NUM
	FILE_TO_ANALYZE="E"$CURR_EMBRYO_NUM
	NUM_FILES_PER_READ=1
	
	# loop over files in each sample (just 1 file / sample for MiSeq, more for HiSeq)
	for ((j=1; j<=$NUM_FILES_PER_READ; j++))
	do

		echo ""
		echo "file in: ./RawData/"$FILE_TO_ANALYZE"_R1_00"$j".fastq"
		echo "file out: ./"$FILE_TO_ANALYZE"_R1_00"$j"_trimmed.fastq" 
		echo ""

		fastq_quality_trimmer -v -t 20 -l 20 -i "./RawData/"$FILE_TO_ANALYZE"_R1_00"$j".fastq" -o "./"$FILE_TO_ANALYZE"_R1_00"$j"_trimmed.fastq" 
		fastq_quality_trimmer -v -t 20 -l 20 -i "./RawData/"$FILE_TO_ANALYZE"_R2_00"$j".fastq" -o "./"$FILE_TO_ANALYZE"_R2_00"$j"_trimmed.fastq" 

	done

	################
	# Part 2: clip #
	################

	echo ""
	echo "now clipping sample : E"$CURR_EMBRYO_NUM

	for ((j=1; j<=$NUM_FILES_PER_READ; j++))
	do
		
		echo ""
		echo "file in: ./"$FILE_TO_ANALYZE"_R1_00"$j"_trimmed.fastq"
		echo "file out: ./ClippedTrimmed/"$FILE_TO_ANALYZE"_R1_00"$j"_cliptrim.fastq"
		echo ""

		fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i "./"$FILE_TO_ANALYZE"_R1_00"$j"_trimmed.fastq" -o "./ClippedTrimmed/"$FILE_TO_ANALYZE"_R1_00"$j"_cliptrim.fastq"
		fastx_clipper -a $ADAPTER_SEQUENCE -l 20 -n -v -i "./"$FILE_TO_ANALYZE"_R2_00"$j"_trimmed.fastq" -o "./ClippedTrimmed/"$FILE_TO_ANALYZE"_R2_00"$j"_cliptrim.fastq"

	done

done

rm *_trimmed.fastq
















