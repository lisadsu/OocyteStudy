#!/bin/bash

# Put .bam files into CuffDiff

SAMPLE_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/HumanHiSeq/"
PICARD_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/picard-tools-1.107/picard-tools-1.107/"
RGLB="libraryName"
RGPL="Illumina"
RGPU="barcode"

declare -a sampleNum=(2)
# 4 5 6 7 9 12 16 17 18 19 20)
cd $SAMPLE_DIR

for i in ${sampleNum[@]}
do

    # Step 2) Label sample with read groups, sample names, etc.

    echo ""
    echo "now executing AddOrReplaceReadGroups on S"$i
    echo ""

    java -Xmx4g -jar $PICARD_DIR"AddOrReplaceReadGroups.jar" INPUT="./Aligned/Bam/S"$i".bam" OUTPUT="./S"$i"_bam/S"$i"_sorted_RG.bam" RGLB=$RGLB RGPL=$RGPL RGPU=$RGPU RGSM=$BASENAME RGID=1 VALIDATION_STRINGENCY=LENIENT

    # Step 3) Reorder reads according to order of reference genome

    echo ""
    echo "now executing ReorderSam on S"$i
    echo ""

    java -Xmx4g -jar $PICARD_DIR"ReorderSam.jar" INPUT="./S"$i"_bam/S"$i"_sorted_RG.bam" OUTPUT="./S"$i"_bam/S"$i"_reorder.bam" REFERENCE=$HUMAN_GENOME_FASTA

    # Step 4) Coordinate sort all reads and create index for .bam file

    echo ""
    echo "now executing SortSam on S"$i
    echo ""

    java -Xmx4g -jar $PICARD_DIR"SortSam.jar" INPUT="./S"$i"_bam/S"$i"_reorder.bam" OUTPUT="./S"$i"_bam/S"$i"_final.bam" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

done

#cuffdiff -o ./CuffDiffOut_v2.0.2_goodPooledOnly -L "good","bad1","bad2","bad3","bad4","bad5","bad6" -b $HUMAN_GENOME_FASTA -u $REF_GTF_HUMAN $SAMPLE_DIR"S2_bam/S2_sorted.bam",$SAMPLE_DIR"S6_bam/S6_sorted.bam",$SAMPLE_DIR"S9_bam/S9_sorted.bam",$SAMPLE_DIR"S17_bam/S17_sorted.bam",$SAMPLE_DIR"S18_bam/S18_sorted.bam",$SAMPLE_DIR"S19_bam/S19_sorted.bam" $SAMPLE_DIR"S4_bam/S4_sorted.bam" $SAMPLE_DIR"S5_bam/S5_sorted.bam" $SAMPLE_DIR"S7_bam/S7_sorted.bam" $SAMPLE_DIR"S12_bam/S12_sorted.bam" $SAMPLE_DIR"S16_bam/S16_sorted.bam" $SAMPLE_DIR"S20_bam/S20_sorted.bam"



