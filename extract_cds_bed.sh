#!/bin/bash

gff=$1
bed=$2

agat_convert_sp_gff2bed.pl --gff $gff --sub CDS --outfile $bed
