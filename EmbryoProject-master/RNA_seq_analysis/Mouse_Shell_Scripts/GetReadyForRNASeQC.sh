#!/bin/bash

# This is just a test script. Use RunRNASeQC.sh to actually go from .bam files to QC report


# Take .bam files that TopHat outputs and run them through RNASeQC
# Define directory with Picard tools

PICARD_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/picard-tools-1.107/picard-tools-1.107/"
RNA_SEQC_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/"
SAMPLE_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/CleanedUp/"

RGLB="libraryName"
RGPL="Illumina"
RGPU="GTGAAA"
RGSM="bad1"
BASENAME="bad1"

# 1) Merge .bam alignment: combine .bam files with aligned and unaligned reads that TopHat outputs
# This one doesn't work for some reason, use samtools merge instead

# java -Xmx4g -jar $PICARD_DIR"AddOrReplaceReadGroups.jar" INPUT="unmapped.bam" OUTPUT="unmapped_RG.bam" RGLB=$RGLB RGPL=$RGPL RGPU=$RGPU RGSM=$RGSM RGID=2 VALIDATION_STRINGENCY=SILENT

# java -Xmx4g -jar $PICARD_DIR"MergeBamAlignment.jar" REFERENCE_SEQUENCE=$MOUSE_GENOME_FASTA ALIGNED_BAM=accepted_hits_RG.bam UNMAPPED_BAM=unmapped_RG.bam OUTPUT=merged_RG.bam PAIRED_RUN=true 

# samtools merge merged.bam accepted_hits.bam unmapped_fixup.bam
 
##################################################
samtools merge $BASENAME"_mergedhits.bam" $SAMPLE_DIR"TophatOut_"$BASENAME"_paired/accepted_hits.bam" $SAMPLE_DIR"TophatOut_"$BASENAME"_singles/accepted_hits.bam"

java -Xmx4g -jar $PICARD_DIR"AddOrReplaceReadGroups.jar" INPUT=$BASENAME"_mergedhits.bam" OUTPUT=$BASENAME"_mergedhits_RG.bam" RGLB=$RGLB RGPL=$RGPL RGPU=$RGPU RGSM=$RGSM RGID=1 
##################################################

# 2) Add read groups

# java -Xmx4g -jar $PICARD_DIR"AddOrReplaceReadGroups.jar" INPUT="merged.bam" OUTPUT="mergedR.bam" RGLB=$RGLB RGPL=$RGPL RGPU=$RGPU RGSM=$RGSM VALIDATION_STRINGENCY=LENIENT

# 3) Sort SAM

# java -Xmx4g -jar $PICARD_DIR"SortSam.jar" INPUT="mergedR.bam" OUTPUT="mergedS.bam" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT

# 4) Reorder SAM since order of contigs in .bam file does not match the order in the .dict file

# FA_IND="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/Mus_musculus_NCBI/Sequence/WholeGenomeFasta/genome.fa.fai"

# java -Xmx4g -jar $PICARD_DIR"ReorderSam.jar" INPUT="mergedS.bam" OUTPUT="mergedReorder.bam" REFERENCE=$MOUSE_GENOME_FASTA VALIDATION_STRINGENCY=SILENT

##################################################
java -Xmx4g -jar $PICARD_DIR"ReorderSam.jar" INPUT=$BASENAME"_mergedhits_RG.bam" OUTPUT=$BASENAME"_mergedhits_reorder.bam" REFERENCE=$MOUSE_GENOME_FASTA
##################################################

# 5) Mark duplicates

# java -Xmx4g -jar $PICARD_DIR"MarkDuplicates.jar" INPUT="mergedReorder.bam" OUTPUT="mergedDuplicates.bam" METRICS_FILE="metricsFile.txt" VALIDATION_STRINGENCY=SILENT

# 6) Index .bam files

# java -Xmx4g -jar $PICARD_DIR"SortSam.jar" INPUT="mergedDuplicates.bam" OUTPUT="mergedFinal.bam" VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate  CREATE_INDEX=true

##################################################
java -Xmx4g -jar $PICARD_DIR"SortSam.jar" INPUT=$BASENAME"_mergedhits_reorder.bam" OUTPUT=$BASENAME"_mergedhits_final.bam" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
##################################################

##################################################
java -Xmx4g -jar $RNA_SEQC_DIR"RNA-SeQC.jar" -r $MOUSE_GENOME_FASTA -t $REF_GTF -o RNA_SeQC_Out -s "1|"$BASENAME"_mergedhits_final.bam|none"
##################################################















