# mk_tests

Scripts to implement an MK test from a set of genome sequence, annotations, 
and variants.

## General overview

1. Obtain the annotation (CDS) coordinates from the GFF using `extract_cds_bed.sh`.
2. Extract consensus for all the coding sequences from all individuals in the VCF/BCF (one consensus 
per haplotype), using `extract_haps_cds.sh`. TODO: for now, it only works on SNP variants.
3. Run `orthofinder` using all the consensus CDSs.
4. TODO

## Python script to extract the CDSs

The program `extract_hap_cds.py` extracts the per-sample coding sequences using 
`samtools faidx` and `bcftools consensus`. The commands process each individual 
CDS entry in the GFF individually, and merges them into the corresponding 
transcripts when necessary.

```
$ python3 extract_hap_cds.py -h
  usage: extract_hap_cds.py [-h] -g GENOME -f GFF -v VCF [-o OUT_DIR] 
                            [-t THREADS] [--snps-only]

  options:
    -h, --help            show this help message and exit
    -g, --genome GENOME   (str) Path to genome in FASTA format.
    -f, --gff GFF         (str) Path to the annotation in FASTA format.
    -v, --vcf VCF         (str) Path to variants in VCF/BCF format.
    -o, --out-dir OUT_DIR
                          (str) Path to output directory [default=.].
    -t, --threads THREADS
                          (int) Number of threads to run in parallel sections 
                          of code [default=1].
    --snps-only           Filter the input variants to only keep SNPs.
```

## Authors

Angel G. Rivera-Colon  
Institute of Ecology and Evolution  
University of Oregon
