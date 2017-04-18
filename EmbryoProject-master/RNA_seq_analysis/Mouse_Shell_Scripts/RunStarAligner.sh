#!/bin/bash
# Type sudo -s before running this script (run as root since memory clearing operation at end requires it)

BASE_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/MouseHiSeq"
STAR_PATH="/home/livia/Desktop/RNA_seq_analysis/STAR_2.3.0e.Linux_x86_64/STAR"
GENOME_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Mus_musculus_NCBI/Sequence/STAR"

# remove all previous STAR directories and make fresh ones

# 1. Generate Genome (only needs to be done once)

#cd $BASE_DIR
#rm -rf "StarOut"
#rm -rf $GENOME_DIR
#mkdir "StarOut"
#mkdir $GENOME_DIR
#cd "StarOut"

#$STAR_PATH --runMode genomeGenerate --genomeDir $GENOME_DIR --genomeFastaFiles $MOUSE_GENOME_FASTA --sjdbGTFfile $REF_GTF_MOUSE --genomeChrBinNbits 12 --runThreadN 2 


#declare -a sampleName=("bad1" "bad2" "good1" "good2" "mixed1" "mixed2")
declare -a sampleName=("bad1" "bad2" "good2" "mixed1" "mixed2")
#declare -a sampleName=("good1")

# loop over samples
for FILENAME in ${sampleName[@]}
do

	echo "Now aligning sample : "$FILENAME	

	cd $BASE_DIR"/StarOut"
	rm -rf "./"$FILENAME
	mkdir "./"$FILENAME 
	cd "./"$FILENAME

	# 2. Align samples

	$STAR_PATH --runMode alignReads --genomeDir $GENOME_DIR --readFilesIn $BASE_DIR"/ClippedTrimmed/"$FILENAME"_R1_stillpaired.fastq" $BASE_DIR"/ClippedTrimmed/"$FILENAME"_R2_stillpaired.fastq" --runThreadN 4 --genomeChrBinNbits 10 --genomeLoad NoSharedMemory --outFilterIntronMotifs RemoveNoncanonicalUnannotated --limitGenomeGenerateRAM 27000000000 

	# 3. Convert SAM to BAM and sort
	cd $BASE_DIR
	mkdir -p "./HTSeq_Count_Out"

	echo "now sorting .bam file and running HTSeq-count"
	
	# convert sam to bam
	samtools view -b -S "./StarOut/"$FILENAME"/Aligned.out.sam" > "./"$FILENAME"_bam/"$FILENAME"_pairedonly.bam"

	# sort bam file
	samtools sort -n "./"$FILENAME"_bam/"$FILENAME"_pairedonly.bam" "./"$FILENAME"_bam/"$FILENAME"_pairedonly_sorted"
	
	# then run HTSeq-count
	python -m HTSeq.scripts.count -f bam "./"$FILENAME"_bam/"$FILENAME"_pairedonly_sorted.bam" $REF_GTF_MOUSE > "./HTSeq_Count_Out/"$FILENAME"_paired_count.txt"

	# clear inactive memory and caches
	sync && echo 3 | tee /proc/sys/vm/drop_caches 

done




