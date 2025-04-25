# mk_tests
Scripts to implement an MK test from a set of genome sequence, annotations, and variants.

## General overview

1. Obtain the annotation (CDS) coordinates from the GFF using `extract_cds_bed.sh`.
2. Extract consensus for all the coding sequences from all individuals in the VCF/BCF (one consensus 
per haplotype), using `extract_haps_cds.sh`. TODO: for now, it only works on SNP variants.
3. Run `orthofinder` using all the consensus CDSs.
4. TODO
