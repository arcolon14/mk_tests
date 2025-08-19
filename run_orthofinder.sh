#!/bin/bash

# Put bash in strict mode
set -e
set -o pipefail

THR=16
NPR=$(($THR/4))

work=/path/to/mk_test
in_data=$work/cds_fasta
of_runs=$work/orthofinder

# List of species for the run
spps=(
    ingroup_1_1
    ingroup_1_2
    ingroup_2_1
    ingroup_2_2
    ingroup_3_1
    ingroup_3_2
    ingroup_4_1
    ingroup_4_2
    ingroup_5_1
    ingroup_5_2
    ingroup_6_1
    ingroup_6_2
    outgroup_1_1
)

# Prepare output directory
N=$(echo ${#spps[@]})
outd=$(date +${of_runs}/%Y%m%d.orthofinder_mktest.N${N})
mkdir -p $outd
cd $outd

# Link the specific input FASTA sequences
cds_dir=$outd/cds
mkdir -p $cds_dir
orgs=$outd/orgs.txt
> $orgs
for sp in ${spps[@]}; do
    # Link the CDS
    ln -s ${in_data}/${sp}.CDS.fa ${cds_dir}/${sp}.fa
    echo $sp >> $orgs
done

# Run OrthoFinder
cmd=(
    orthofinder
    -t $THR
    -a $NPR
    -M "msa"
    -d # Input is DNA
    -A "mafft"
    -T "fasttree"
    -o $outd/mktest.N${N}
    -y
    -f $cds_dir
)
echo "${cmd[@]}"
"${cmd[@]}"
