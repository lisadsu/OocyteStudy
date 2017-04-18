#!/bin/bash

# Run RNASeQC and also prepare files for differential expression analysis
# Works only on mapped reads -- does not take into account unmapped reads. Stats for % of reads mapped can be found in tophat output directories

RNA_SEQC_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/"
PICARD_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/picard-tools-1.107/picard-tools-1.107/"
SAMPLE_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/CleanedUp/"

RGLB="libraryName"
RGPL="Illumina"
RGPU="barcode"

declare -a arr=("bad1" "bad2" "good1" "good2" "mixed1" "mixed2" "pos1" "pos2")

for BASENAME in ${arr[@]}
do

    # Step 1) Take .bam files that TopHat outputs and merge mapped paired end and single reads (not unmapped reads since there is some problem with the formatting

    samtools merge $BASENAME"_mergedhits.bam" $SAMPLE_DIR"TophatOut_"$BASENAME"_paired/accepted_hits.bam" $SAMPLE_DIR"TophatOut_"$BASENAME"_singles/accepted_hits.bam"

    # Step 2) Label sample with read groups, sample names, etc.

    java -Xmx4g -jar $PICARD_DIR"AddOrReplaceReadGroups.jar" INPUT=$BASENAME"_mergedhits.bam" OUTPUT=$BASENAME"_mergedhits_RG.bam" RGLB=$RGLB RGPL=$RGPL RGPU=$RGPU RGSM=$BASENAME RGID=1 

    # Step 3) Reorder reads according to order of reference genome

    java -Xmx4g -jar $PICARD_DIR"ReorderSam.jar" INPUT=$BASENAME"_mergedhits_RG.bam" OUTPUT=$BASENAME"_mergedhits_reorder.bam" REFERENCE=$MOUSE_GENOME_FASTA

    # Step 4) Coordinate sort all reads and create index for .bam file

    java -Xmx4g -jar $PICARD_DIR"SortSam.jar" INPUT=$BASENAME"_mergedhits_reorder.bam" OUTPUT=$BASENAME"_mergedhits_final.bam" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

    # Step 5) Run RNA-SeQC and generate report

    java -Xmx4g -jar $RNA_SEQC_DIR"RNA-SeQC.jar" -r $MOUSE_GENOME_FASTA -t $REF_GTF -o "RNA_SeQC_Out/"$BASENAME -s "1|"$BASENAME"_mergedhits_final.bam|none"

    mkdir -p $BASENAME"_bam/"
    mv *_final.bam $BASENAME"_bam/"
    mv *_final.bai $BASENAME"_bam/"
    rm *_RG.bam
    rm *_reorder.bam
    rm *_mergedhits.bam

done



