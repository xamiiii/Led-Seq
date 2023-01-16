#!/bin/bash

library=$1
ANNOTATION=$2


# overlap if 1 bp of read overlaps with transcript annotation
bedtools intersect -loj -s -a "${library}.sorted.dedup.primary.bed" -b $ANNOTATION > "${library}.sorted.dedup.primary.intersect.loj"
awk 'BEGIN{lastread=""} {if ($4==lastread) next; print $0; lastread=$4}' "${library}.sorted.dedup.primary.intersect.loj" > "${library}.sorted.dedup.primary.intersect.loj.singular"   #always take the first intersection if read has multiple hits in annotations

rm "${library}.sorted.dedup.primary.intersect.loj"

TOTAL_READS=$(wc -l "${library}.sorted.dedup.primary.bed"| awk '{print $1}')
echo "total reads: " ${TOTAL_READS} > "${library}.composition"
PROTEIN_CODING=$(grep gene_biotype=protein_coding "${library}.sorted.dedup.primary.intersect.loj.singular" | wc -l)
echo "protein_coding total: " ${PROTEIN_CODING} >> "${library}.composition"
PROTEIN_CODING_PERCENT=$(python3 -c "print(round(${PROTEIN_CODING}/${TOTAL_READS}*100,1))")
echo "protein_coding [%]  : " ${PROTEIN_CODING_PERCENT} >> "${library}.composition"

NCRNA=$(grep -E 'gene_biotype=ncRNA|RNase_P_RNA|tmRNA|gene_biotype=antisense_RNA|SRP_RNA|ssrS' "${library}.sorted.dedup.primary.intersect.loj.singular" | wc -l)
echo "ncRNA total: " ${NCRNA} >> "${library}.composition"
NCRNA_PERCENT=$(python3 -c "print(round(${NCRNA}/${TOTAL_READS}*100,1))")
echo "ncRNA [%]  : " ${NCRNA_PERCENT} >> "${library}.composition"

RRNA=$(grep -E '5S_ribosomal_RNA|23S_ribosomal_RNA|16S_ribosomal_RNA|gene_biotype=rRNA' "${library}.sorted.dedup.primary.intersect.loj.singular" | wc -l)
echo "rRNA total: " ${RRNA} >> "${library}.composition"
RRNA_PERCENT=$(python3 -c "print(round(${RRNA}/${TOTAL_READS}*100,1))")
echo "rRNA [%]  : " ${RRNA_PERCENT} >> "${library}.composition"

TRNA=$(grep -E 'tRNA' "${library}.sorted.dedup.primary.intersect.loj.singular" | wc -l)
echo "tRNA total: " ${TRNA} >> "${library}.composition"
TRNA_PERCENT=$(python3 -c "print(round(${TRNA}/${TOTAL_READS}*100,1))")
echo "tRNA [%]  : " ${TRNA_PERCENT} >> "${library}.composition"

NOT_ANNOTATED=$(awk '$10=="-1"' "${library}.sorted.dedup.primary.intersect.loj.singular" | wc -l)
echo "not_annotated total: " ${NOT_ANNOTATED} >> "${library}.composition"
NOT_ANNOTATED_PERCENT=$(python3 -c "print(round(${NOT_ANNOTATED}/${TOTAL_READS}*100,1))")
echo "not_annotated [%]  : " ${NOT_ANNOTATED_PERCENT} >> "${library}.composition"

rm "${library}.sorted.dedup.primary.intersect.loj.singular"
