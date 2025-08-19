#!/bin/bash
# Put bash in strict mode
set -e
set -o pipefail

THR=12
# This working directory contains the outputs of Orthofinder
work=/path/to/mk_test
# Single copy orthogroup sequences from Orthofinder
indir=$work/orthofinder_out/Single_Copy_Orthologue_Sequences/
# List of single copy orthogroup ids from Orthofinder
sco_list=$work/orthofinder_out/Orthogroups/Orthogroups_SingleCopyOrthologues.txt
outdir=$work/out_refined_msa

cat $sco_list | grep -v '^#' |
while read og_id; do
    # Align the single copy ortholog sequences with macse
    in_fasta=$indir/${og_id}.fa
    # Outputs
    out_aa=$outdir/${og_id}_AA.fa
    out_nt=$outdir/${og_id}_NT.fa
    log=$outdir/${og_id}.macse.log
    cmd=(
        macse
        -prog alignSequences
        -seq $in_fasta
        -out_NT $out_nt
        -out_AA $out_aa
    )
    echo "${cmd[@]} > ${log}"
done | parallel -j $THR

# Loop back over the resulting alignments and extract the outgroup
# TODO: Probably cleaner in a python script
splitdir=$work/split_msa
# ID of orthogroup as it appears on FASTA
outgroup_id="outgroup"
cat $sco_list | grep -v '^#' |
while read og_id; do
    in_msa=$outdir/${og_id}_NT.fa
    # Extract ingroup
    ingrp_msa=$splitdir/${og_id}_ingroup.fa
    echo "cat $in_msa | seqtk seq -l0 | paste - - | grep -v ${outgroup_id} | tr '\t' '\n' | seqtk seq -l60 > $ingrp_msa"
    # Extract outgroup
    outgrp_msa=$splitdir/${og_id}_outgroup.fa
    echo "cat $in_msa | seqtk seq -l0 | paste - - | grep ${outgroup_id} | tr '\t' '\n' | seqtk seq -l60 > $outgrp_msa"
done | parallel -j $THR
