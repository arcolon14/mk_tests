#!/usr/bin/env python3
import sys, os, gzip, argparse, datetime, re, subprocess
from shutil import which
from multiprocessing import Pool

# Some constants
PROG = sys.argv[0].split('/')[-1]
SNPS_ONLY = False
FA_LINE_WIDTH = 60

def parse_args(prog=PROG):
    '''Set and verify command line options.'''
    p = argparse.ArgumentParser()
    p.add_argument('-g', '--genome', required=True, 
                   help='(str) Path to genome in FASTA format.')
    p.add_argument('-f', '--gff', required=True, 
                   help='(str) Path to the annotation in GFF format.')
    p.add_argument('-v', '--vcf', required=True, 
                   help='(str) Path to variants in VCF/BCF format.')
    p.add_argument('-o', '--out-dir', required=False, default='.',
                   help='(str) Path to output directory [default=.].')
    p.add_argument('-t', '--threads', required=False, type=int, default=1,
                   help='(int) Number of threads to run in parallel sections of code [default=1].')
    p.add_argument('--snps-only', action='store_true',
                   help='Filter the input variants to only keep SNPs.')
    # Check inputs
    args = p.parse_args()
    assert os.path.exists(args.genome)
    assert os.path.exists(args.gff)
    assert os.path.exists(args.vcf)
    args.out_dir = args.out_dir.rstrip('/')
    assert args.threads >= 1
    # Set some constants
    SNPS_ONLY = args.snps_only
    return args


def date() -> str:
    '''Print the current date in YYYY-MM-DD format.'''
    return datetime.datetime.now().strftime("%Y-%m-%d")

def time() -> str:
    '''Print the current time in HH:MM:SS format.'''
    return datetime.datetime.now().strftime("%H:%M:%S")

def check_executables() -> None:
    '''Check that the required executables are in the PATH.'''
    print('\nChecking for required executables...', flush=True)
    for program in ['bcftools', 'samtools']:
        if which(program) is None:
            sys.exit(f"Error: \'{program}\' not found in PATH.")
        else:
            print(f'    {program} found in PATH.')

def check_input_fasta(fasta_f: str) -> None:
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
        subprocess.run(['samtools', 'faidx', fasta_f], check=True)

def check_input_vcf(vcf_f: str) -> None:
    '''Check input VCF/BCF and confirm that it is indexed.
    Args:
        vcf_f : (str) Path to varints file in VCF or BCF format.
    Returns:
        None
    '''
    print('\nChecking input VCF/BCF and index...', flush=True)
    if os.path.exists(f'{vcf_f}.csi'):
        print('    VCF/BCF CSI index found.')
    else:
        print('    VCF/BCF CSI index not found. Generaring it with `bcftools index`.')
        subprocess.run(['bcftools', 'index', vcf_f], check=True)

def extract_vcf_samples(vcf_f: str) -> list:
    '''
    Extract the samples from the VCF header.
    Args:
        vcf_f : (str) Path to input VCF/BCF
    Returns:
        samples : list of samples
    '''
    print('\nExtracting samples from the VCF/BCF header...', flush=True)
    cmd = ['bcftools', 'query', '--list-samples', vcf_f]
    process = subprocess.run(cmd, capture_output=True, check=True, text=True)
    stdout = process.stdout
    samples = stdout.strip('\n').split('\n')
    print(f'    Found {len(samples):,} samples from the input VCF/BCF.', flush=True)
    return samples

class Transcript:
    '''
    Class to store the mRNA/transcript records from the GFF file.
    '''
    def __init__(self, id: str, chrom: str, start: int, stop: int, strand: str):
        assert type(start) == int
        assert type(stop) == int
        assert stop >= start
        assert strand in {'-', '+'}
        # Attributes
        self.id = id
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.strand = strand
        # Set empty list for the exons
        self.exons_l = list()
    def __str__(self):
        return f'{self.id} {self.chrom} {self.start} {self.stop} {self.strand}'

class CodingExon:
    '''
    Class to store the CDS entries from the GFF file.
    '''
    def __init__(self, id: str, chrom: str, start: int, stop: int, phase: int):
        assert type(start) == int
        assert type(stop) == int
        assert stop >= start
        assert type(phase) == int
        assert 0 <= phase < 3 # Phases can only be 0,1,2
        # Attributes
        self.id = id
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.phase = phase
    def __str__(self):
        return f'{self.id} {self.chrom} {self.start} {self.stop} {self.phase}'

def load_cds_from_gff(gff_f: str) -> dict:
    '''
    Parse the GFF and extract the coding sequences.
    Args:
        gff_f : (str) Path to annotations in GFF format.
    Returns:
        annotations : (dict) Dictionary of Transcript objects.
    '''
    print('\nParsing the GFF and extracting annotations...', flush=True)
    # Main output
    annotations = dict()
    # Parse the GFF
    records = 0
    n_transc = 0
    n_cds = 0
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
                # If the element is an mRNA and trasncript, then this is 
                # the first layer of the annotation.
                id = None
                attributes = fields[8].split(';')
                for attribute in attributes:
                    if attribute.startswith('ID'):
                        # Extract the ID
                        # This is done per-attribute because some of the 
                        # name fields have white spaces and other characters.
                        pairs = re.split('[ =]', attribute)
                        assert len(pairs) == 2, f'Error: attributes not being parsed correctly as key=value pairs:\n{pairs}\n{fields[8]}'
                        # The ID value is the second element of the pair
                        id = pairs[1].rstrip(' \"\'').lstrip(' \"\'')
                # Prepare the Transcript class
                transcript = Transcript(id, fields[0], int(fields[3]), 
                                        int(fields[4]), fields[6])
                # Store in the dictionary
                if id in annotations:
                    sys.exit(f"Error: Transcript id {id} is not unique.")
                annotations[id] = transcript
                n_transc += 1
            elif fields[2] == 'CDS':
                # Process the CDSs. These are the second layer of the 
                # annotation and must be added to the corresponding 
                # transcript.
                id = None
                parent = None
                attributes = fields[8].split(';')
                for attribute in attributes:
                    if attribute.startswith('ID'):
                        # Extract the ID
                        # Again, this is done per-attribute because some of 
                        # the name fields have white spaces and other characters.
                        pairs = re.split('[ =]', attribute)
                        assert len(pairs) == 2, f'Error: attributes not being parsed correctly as key=value pairs:\n{pairs}\n{fields[8]}'
                        # The ID value is the second element of the pair
                        id = pairs[1].rstrip(' \"\'').lstrip(' \"\'')
                    elif attribute.startswith('Parent'):
                        # Extract the Parent
                        pairs = re.split('[ =]', attribute)
                        assert len(pairs) == 2, f'Error: attributes not being parsed correctly as key=value pairs:\n{pairs}\n{fields[8]}'
                        # The ID value is the second element of the pair
                        parent = pairs[1].rstrip(' \"\'').lstrip(' \"\'')
                # Generate the CodingExon object
                cds = CodingExon(id, fields[0], int(fields[3]), 
                                        int(fields[4]), int(fields[7]))
                # Add this to the existing Transcript object
                assert parent in annotations, f'Error: CDS {id} seen before parent transcript {parent}.'
                assert isinstance(annotations[parent], Transcript)
                annotations[parent].exons_l.append(cds)
                n_cds += 1
            else:
                continue
    # Report to log
    print(f'    Read {records:,} records from the GFF file.')
    print(f'    Extracted {n_transc:,} transcripts composed of {n_cds:,} coding sequences.', flush=True)
    return annotations

def rev_comp(sequence: str) -> str:
    '''
    Reverser compliment a nucleotide sequence
    Args:
        sequence: (str) input nucleotide sequence
    Returns:
        rev_seq: (str) reverse complimented sequence
    '''
    rev = list()
    for nt in sequence.upper():
        if nt == 'A':
            rev.append('T')
        elif nt == 'C':
            rev.append('G')
        elif nt == 'G':
            rev.append('C')
        elif nt == 'T':
            rev.append('A')
        else:
            rev.append('N')
    rev_seq = ''.join(rev[::-1])
    return rev_seq

def process_sample(sample: str, haplotype: int, annotations: dict,
                   genome_f: str, vcf_f: str, out_dir: str='.',
                   fa_line_width=FA_LINE_WIDTH) -> None:
    '''
    Process the data for a single sample using the extracted 
    annotations and variant file.
    Args:
        sample: (str) current sample to process
        haplotype: (int) which haplotype to extract variants for
        annotations: (dict) dictionary of transcript and cds objects
        genome_f: (str) Path to genome in FASTA format
        vcf_f: (str) Path to variants in VCF/BCF format
        out_dir: (str) Path to directory where to store outputs [default=.]
    Returns:
        None
    '''
    print(f'    Working on {sample}, haplotype {haplotype}')
    # Prepare the output
    ouf_fa = f'{out_dir}/{sample}_{haplotype}.CDS.fa'
    with open(ouf_fa, 'w') as fh:
        # Loop over the annotations
        for j, trans_id in enumerate(annotations):
            # Select the target transcript
            transcript = annotations[trans_id]
            # Process that transcript
            cds_seq = process_transcript(sample, haplotype, transcript,
                                         genome_f, vcf_f)
            # Save to the file
            fh.write(f'>{trans_id}\n')
            # Wrap the sequence lines up to `fa_line_width` characters
            for start in range(0, len(cds_seq), fa_line_width):
                seq_line = cds_seq[start:(start+fa_line_width)]
                fh.write(f'{seq_line}\n')

def process_transcript(sample: str, haplotype: int, transcript: Transcript, 
                       genome_f: str, vcf_f: str) -> str:
    '''
    Process the the coding sequences and variants for a single transcript.
    Args:
        sample: (str) current sample to process
        haplotype: (int) which haplotype to extract variants for
        transcript: (Transcript) a transcript object
        genome_f: (str) Path to genome in FASTA format
        vcf_f: (str) Path to variants in VCF/BCF format
    Returns:
        cds_seq: (str) Extracted CDS sequence with individual variants
    '''
    assert isinstance(transcript, Transcript)
    assert haplotype in {1,2}
    cds_seq = ''
    # Loop over the individual CDSs in the transcript
    for cds in transcript.exons_l:
        assert isinstance(cds, CodingExon)
        # Extract a CDS with the variants for that haplotype
        var_seq = extract_cds(sample, haplotype, cds, genome_f, vcf_f)
        # Append to the main sequence
        cds_seq += var_seq
    # Process if reverse complimented
    if transcript.strand == '-':
        cds_seq = rev_comp(cds_seq)
    return cds_seq

def extract_cds(sample: str, haplotype: int, cds: CodingExon, 
                genome_f: str, vcf_f: str) -> str:
    '''
    Process the the coding sequences and variants for a single transcript.
    Args:
        sample: (str) current sample to process
        haplotype: (int) which haplotype to extract variants for
        cds: (CodingExon) object for the individual CDS
        genome_f: (str) Path to genome in FASTA format
        vcf_f: (str) Path to variants in VCF/BCF format
    Returns:
        var_seq: (str) Extracted sequence with individual variants
    '''
    assert isinstance(cds, CodingExon)
    var_seq = ''
    # 1. Prepare and run the samtools command to extract the reference
    samt_cmd = ['samtools', 'faidx', 
                f'{genome_f}', f'{cds.chrom}:{cds.start}-{cds.stop}']
    samt_proc = subprocess.Popen(samt_cmd,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 text=True)
    # Check output
    # This sometimes returns None (not 0) on success, so adding the
    # additional check to be sure.
    if samt_proc.returncode is not None and samt_proc.returncode != 0:
        print(samt_proc.returncode)
        samt_str = ' '.join(samt_cmd)
        sys.exit(f"Error: `{samt_str}` exited with non-zero status.")

    # 2. Prepare and run the BCFtools command to add the variants on 
    # the extracted sequence.
    # Note: The `regions-overlap 0` is added to prevent issues with 
    # indels overlapping the boundaries of the target sequences, as
    # defined by the BCFtools documentation.
    bcft_cmd = ['bcftools', 'consensus', 
                '--missing', 'N',
                '--regions-overlap', '0',
                '--samples', f'{sample}',
                '--haplotype', f'{haplotype}']
    # Add conditional arguments
    if SNPS_ONLY:
        bcft_cmd.append(['--exclude', 'TYPE=\'snp\''])
    bcft_cmd.append(f'{vcf_f}')
    bcft_proc = subprocess.Popen(bcft_cmd,
                                 stdin=samt_proc.stdout,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 text=True)
    # Process the stdout
    output, errors = bcft_proc.communicate()
    # Check output
    if bcft_proc.returncode != 0:
        cmd_str = ' '.join(samt_cmd) + ' | ' + ' '.join(bcft_cmd)
        sys.exit(f"Error: `{cmd_str}` exited with non-zero status.\n{errors}")
    # Extract the stdout as a sequence
    stdout_l = output.strip('\n').split('\n')
    for element in stdout_l:
        if not element.startswith('>'):
            var_seq += element
    # Make sure that sequence of non-zero length is returned
    if len(var_seq) < 1:
        sys.exit(f"Error: {sample}_{haplotype} consensus for {cds.chrom}:{cds.start}-{cds.stop} returned a sequence of length 0.")
    var_seq = var_seq.upper()
    return var_seq

def process_all_samples(samples: list, annotations: dict, 
                        genome_f: str, vcf_f: str, out_dir: str='.', 
                        threads:int=1, ploidy:int=2) -> None:
    '''
    Process the data for all the samples using the extracted annotations
    and variant file.
    Args:
        samples: (list) list of samples
        annotations: (dict) dictionary of transcript and cds objects
        genome_f: (str) Path to genome in FASTA format
        vcf_f: (str) Path to variants in VCF/BCF format
        out_dir: (str) Path to directory where to store outputs [default=.]
        threads: (str) Number of threads to run samples in parallel [default=1]
        ploidy: (str) How many haplotypes to extract per sample [default=2]
    Returns:
        None
    '''
    print(f'\nProcessing samples...', flush=True)
    # When working with a single thread...
    if threads == 1:
        # Iterate over samples one at a time
        for sample in samples:
            for h in range(0, ploidy):
                hap = h+1
                process_sample(sample, hap, annotations, 
                               genome_f, vcf_f, out_dir)
    else:
        # When working multithreaded...
        # First, adjust the number of threads if needed
        if threads > len(samples)*ploidy:
            threads = len(samples)*ploidy
        # Call `process_sample` multithreaded
        with Pool(threads) as pool:
            args = list()
            for sample in samples:
                for h in range(0, ploidy):
                    hap = h+1
                    task_args = (sample, hap, annotations, 
                                 genome_f, vcf_f, out_dir)
                    args.append(task_args)
            pool.starmap(process_sample, args)

def main():
    print(f'{PROG} started on {date()} {time()}.')
    # Parse args
    args = parse_args()
    # Check that the needed programs are available
    check_executables()
    # Check the input fasta
    check_input_fasta(args.genome)
    # Check the input VCF/BCF
    check_input_vcf(args.vcf)
    # Extract the samples from the VCF file
    samples = extract_vcf_samples(args.vcf)
    # Extract the CDS annotations from the GFF
    annotations = load_cds_from_gff(args.gff)
    # Process all samples
    process_all_samples(samples, annotations, args.genome, args.vcf, args.out_dir, args.threads)
    # Done!
    print(f'\n{PROG} finished on {date()} {time()}.')

# Run Code
if __name__ == '__main__':
    main()
