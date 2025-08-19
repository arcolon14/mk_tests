#!/bin/bash

gff=$1
bed=$2

cat $gff | grep -v '^#' | \
    awk '$3 == "CDS" {print $1, $4-1, $5, $9, 0, $7}' | \
    tr ' ' '\t' > $bed

# or...

agat_convert_sp_gff2bed.pl --gff $gff --sub CDS --outfile ${bed}12
