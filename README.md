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

```sh
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

## Calculate the codong positions of each coding site

The program `position_in_codon.py` parses an annotation file (GFF) and for 
each coding site, calculate the the codon number (i.e., which codon in amino 
acid sequence the site belongs to) and codon position (i.e., 1st, 2nd, or 3rd 
position of the codon).

```sh
$ python3 position_in_codon.py -h
  position_in_codon.py started on 2025-05-08 23:47:18.
  usage: position_in_codon.py [-h] -g GFF [-o OUT_DIR]

  Parse an annotation file (GFF) and for each coding site, calculate the the codon
  number (i.e., which codon in amino acid sequence the site belongs to) and codon
  position (i.e., 1st, 2nd, or 3rd position of the codon).
  
  options:
    -h, --help            show this help message and exit
    -g GFF, --gff GFF     (str) Path to the annotation in GFF format.
    -o OUT_DIR, --out-dir OUT_DIR
                          (str) Path to output directory [default=.].
```

The script generates the output table, `position_in_codons.tsv`. 
Output example:

```sh
#Chrom  Bp     Transcript  CDS    CodonNum  PosInCodon
chr01   53103  mrna-1      cds-1  1         1
chr01   53104  mrna-1      cds-1  1         2
chr01   53105  mrna-1      cds-1  1         3
chr01   53106  mrna-1      cds-1  2         1
chr01   53107  mrna-1      cds-1  2         2
chr01   53108  mrna-1      cds-1  2         3
chr01   53109  mrna-1      cds-1  3         1
chr01   53110  mrna-1      cds-1  3         2
chr01   53111  mrna-1      cds-1  3         3
...
```

## Authors

Angel G. Rivera-Colon  
Institute of Ecology and Evolution  
University of Oregon
