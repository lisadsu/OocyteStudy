#!/bin/bash
# Put .bam files into CuffDiff

SAMPLE_BASE_DIR="/host/Users/Livia/Desktop/IVF/RnaSeqAnalysis/RawData/MouseHiSeq/"
cd $SAMPLE_BASE_DIR

# Put .bam files into CuffDiff 2.0.2
# but only sorted paired end ones

cuffdiff -o ./CuffDiffOut_v2.0.2_justgoodbad_concordantpairsonly -L bad,good -p 4 -b $MOUSE_GENOME_FASTA -u $REF_GTF_MOUSE $SAMPLE_BASE_DIR"bad1_bam/bad1_pairedonly.bam",$SAMPLE_BASE_DIR"bad2_bam/bad2_pairedonly.bam" $SAMPLE_BASE_DIR"good1_bam/good1_pairedonly.bam",$SAMPLE_BASE_DIR"good2_bam/good2_pairedonly.bam"

# Put .bam files into CuffDiff 2.1.1
# but only sorted paired end ones

/home/livia/Desktop/RNA_seq_analysis/cufflinks-2.1.1/cuffdiff -o ./CuffDiffOut_v2.1.1_justgoodbad_concordantpairsonly -L bad,good -p 4 -b $MOUSE_GENOME_FASTA -u $REF_GTF_MOUSE $SAMPLE_BASE_DIR"bad1_bam/bad1_pairedonly.bam",$SAMPLE_BASE_DIR"bad2_bam/bad2_pairedonly.bam" $SAMPLE_BASE_DIR"good1_bam/good1_pairedonly.bam",$SAMPLE_BASE_DIR"good2_bam/good2_pairedonly.bam"

