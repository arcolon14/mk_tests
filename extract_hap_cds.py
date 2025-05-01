#!/usr/bin/env python3
import sys, os, gzip, argparse, datetime
from subprocess import run
from shutil import which

# Some constants
PROG = sys.argv[0].split('/')[-1]

def parse_args(prog=PROG):
    '''Set and verify command line options.'''
    p = argparse.ArgumentParser()
    p.add_argument('-g', '--genome', required=True, 
                   help='(str) Path to genome in FASTA format.')
    p.add_argument('-f', '--gff', required=True, 
                   help='(str) Path to the annotation in FASTA format.')
    p.add_argument('-v', '--vcf', required=True, 
                   help='(str) Path to variants in VCF/BCF format.')
    p.add_argument('-o', '--out-dir', required=False, default='.',
                   help='(str) Path to output directory [default=.].')
    p.add_argument('--snps-only', action='store_true',
                   help='Filter the input variants to only keep SNPs.')
    p.add_argument('--biallelic-only', action='store_true',
                   help='Filter the input variants to only keep biallelic sites.')
    # Check inputs
    args = p.parse_args()
    assert os.path.exists(args.genome)
    assert os.path.exists(args.gff)
    assert os.path.exists(args.vcf)
    args.out_dir = args.out_dir.rstrip('/')
    return args


def date():
    '''Print the current date in YYYY-MM-DD format.'''
    return datetime.datetime.now().strftime("%Y-%m-%d")

def time():
    '''Print the current time in HH:MM:SS format.'''
    return datetime.datetime.now().strftime("%H:%M:%S")

def check_executables():
    '''Check that the required executables are in the PATH.'''
    print('\nChecking for required executables...', flush=True)
    for program in ['bcftools', 'samtools']:
        if which(program) is None:
            sys.exit(f"Error: \'{program}\' not found in PATH.")
        else:
            print(f'    {program} found in PATH.')

def check_input_fasta(fasta_f):
    '''Check input FASTA and confirm that it is indexed.
    Args:
        fasta_f : (str) Path to FASTA file.
    Returns:
        None
    '''
    print('\nChecking input FASTA and index...', flush=True)
    if os.path.exists(f'{fasta_f}.fai'):
        print('    FASTA FAI index found.')
    else:
        print('    FASTA FAI index not found. Generaring it with `samtools faidx`.')
        run(['samtools', 'faidx', fasta_f], check=True)

def extract_vcf_samples(vcf_f):
    '''
    Extract the samples from the VCF header.
    Args:
        vcf_f : (str) Path to input VCF/BCF
    Returns:
        samples : list of samples
    '''
    print('\nExtracting samples from the VCF/BCF header...', flush=True)
    cmd = ['bcftools', 'query', '--list-samples', vcf_f]
    process = run(cmd, capture_output=True, check=True, text=True)
    stdout = process.stdout
    samples = stdout.strip('\n').split('\n')
    print(f'    Found {len(samples):,} from the input VCF/BCF.', flush=True)
    return samples

def load_cds_from_gff(gff_f):
    '''
    Parse the GFF and extract the coding sequences.
    Args:
        gff_f : (str) path to annotations in GFF format.
    Return:
        annotations : (dict)
    '''
    print('\nParsing the GFF and extracting annotations...', flush=True)
    # Main output
    annotations = dict()
    # Parse the GFF
    records = 0
    with open(gff_f) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            line = line.strip('\n')
            if len(line) == 0:
                continue
            # Now process the records
            records += 1
            fields = line.split('\t')
            # Process the different elements
            if fields[2] in {'mRNA', 'transcript'}:
                # If the element is an mRNA and trasncript extract the 
                # ID. This is the key for the dictionary.
                attributes = fields[8].split(';')
                for attribute in attributes:
                    pairs = attribute.split('=')
                    key = pairs[0].rstrip(' ').lstrip(' ')
                    value = pairs[1].rstrip(' ').lstrip(' ')
                    # Store the ID as the key of the dictionary
                    if key == 'ID':
                        annotations.setdefault(value, [])
            elif fields[2] == 'CDS':
                print(fields)
            else:
                continue


            if records > 50:
                break
    print(annotations)

def main():
    print(f'{PROG} started on {date()} {time()}.')
    # Parse args
    args = parse_args()
    # Check that the needed programs are available
    check_executables()
    # Check the input fasta
    check_input_fasta(args.genome)
    # Extract the samples from the VCF file
    samples = extract_vcf_samples(args.vcf)
    # Extract the CDS annotations from the GFF
    annotations = load_cds_from_gff(args.gff)
    # Process all samples
    # process_all_samples(samples, args.genome, args.gf)
    # Done!
    print(f'\n{PROG} finished on {date()} {time()}.')

# Run Code
if __name__ == '__main__':
    main()
