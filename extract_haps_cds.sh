#!/bin/bash
set -e
set -o pipefail

work=/sietch_colab/ariverac/balanus_genome/mk_test

# Input sequences
ref=$work/ref_data
genome=$ref/BalGla.fa
# List of scaffolds to process. This is taken from the FAI
# file. For now, are used to access the chromosome-
# specific BCFs.
scaffs=( $(cat ${genome}.fai | cut -f1) )
# Variants dir
# vars=$work/variants/ref_indv
vars=$work/variants/wgs_n6
# BED file containing the annotations spans
bed=$ref/BalGla.CDS.bed12
# Output directory
outp=$work/hap_cds
# List of individuals to process, this could be a
# separate file, but we are taking it from the BCF header.
samples=( $(bcftools view -h $vars/Bgland_${scaffs[0]}.filtered.bcf | \
            grep '^#CHROM' | cut -f 10- | tr '\t' '\n') )

# Clear any temporary directories. This prevents some issues
# with existing FAIs, when rerunning.
if [ -d $outp/tmp ]
then
    rm -r $outp/tmp
fi
# Then create the tmp directory
mkdir $outp/tmp


# First, loop over each individual
for sample in ${samples[@]}
do
    # Then loop over each haplotype
    for hap in {1..2}
    do
        echo "Working on ${sample}, haplotype ${hap}..."
        # Initialize the output genome FASTA for each haplotype
        # This is a temporary file
        out_hap_fa=$outp/tmp/${sample}_${hap}.fa
        > $out_hap_fa
        # Initalize the output CDS FASTA per each haplotype
        out_cds_fa=$outp/${sample}_${hap}.CDS.fa
        > $out_cds_fa
        # Then loop over each chromosome
        for scaff in ${scaffs[@]}
        do
            echo "    Processing ${scaff}"
            # Select the target variants for that chroomosome
            bcf=$vars/Bgland_${scaff}.filtered.bcf

            # Right now, we only have variants for scaffs
            # > 1Mbp, so some BCFs might not exists. Handle
            # accordingly.
            if [ -e $bcf ]
            then
                # If the BCF exists...
                # Extract the target chromosome and add variants.
                # The seqtk seq -U is for converting to uppercase for
                # later compatibility. For now, indels are removed, but
                # we should probably check the --mark-ins and --mark-del
                # options.
                samtools faidx $genome "${scaff}" | \
                    bcftools consensus --samples ${sample} \
                        --exclude "TYPE!='snp'" \
                        --haplotype ${hap} \
                        --missing "N" \
                        $bcf | \
                    seqtk seq -U | \
                    fold -w 60 >> $out_hap_fa
            else
                # If the BCF does not exists, don't extract
                # any variants. Just append reference sequence.
                echo "        Variants for ${scaff} not found. Appending reference."
                samtools faidx $genome "${scaff}" | \
                    seqtk seq -U | \
                    fold -w 60 >> $out_hap_fa
            fi

        done
        # Extract the CDS spans based on the modified sequence for
        # that individual.
        # Note: I don't think this will work with indels.
        bedtools getfasta \
                -fi $out_hap_fa \
                -bed $bed \
                -nameOnly \
                -s \
                -split | \
            fold -w 60 >> $out_cds_fa
        n_cds=$(cat $out_cds_fa | grep -c '^>')
        echo -e "\n    Extracted ${n_cds} sequences.\n"
    # Remove the tmp files for that haplotype
    rm $out_hap_fa*
    done
done
