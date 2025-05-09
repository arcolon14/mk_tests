#!/usr/bin/env python3
import sys, os, datetime, argparse, re
PROG = sys.argv[0].split('/')[-1]

def parse_args(prog=PROG):
    '''Set and verify command line options.'''
    p = argparse.ArgumentParser()
    p.add_argument('-g', '--gff', required=True, 
                   help='(str) Path to the annotation in GFF format.')
    p.add_argument('-o', '--out-dir', required=False, default='.',
                   help='(str) Path to output directory [default=.].')
    # Check inputs
    args = p.parse_args()
    assert os.path.exists(args.gff)
    args.out_dir = args.out_dir.rstrip('/')
    return args

def date() -> str:
    '''Print the current date in YYYY-MM-DD format.'''
    return datetime.datetime.now().strftime("%Y-%m-%d")

def time() -> str:
    '''Print the current time in HH:MM:SS format.'''
    return datetime.datetime.now().strftime("%H:%M:%S")

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

def calculate_codon_positions(annotations:dict, outdir:str='.') -> None:
    '''
    Calculate the per-site position across codons from the annotations.
    Args:
        annotations: (dict) Dictionary with the annotations
        outdir: (str) Path to output directory
    Returns:
        None
    '''
    print('\nCalculating the position in the exon for each site in the CDS...')
    total_sites = 0
    # Prepare output
    with open(f'{outdir}/position_in_codons.tsv', 'w') as fh:
        fh.write('#Chrom\tBp\tTranscript\tCDS\tCodonNum\tPosInCodon\n')
        # Loop over the transcripts and CDSs
        for trans_id in annotations:
            transcript = annotations[trans_id]
            codon_num = 0
            # Handle genes in positive strand
            if transcript.strand == '+':
                for cds in transcript.exons_l:
                    # Loop over the range of the CDS
                    for i, pos in enumerate(range(cds.start-1, cds.stop)):
                        codon_pos = ((i-cds.phase)%3)+1
                        if codon_pos == 1:
                            codon_num+=1
                        row = f'{transcript.chrom}\t{pos+1}\t{transcript.id}\t{cds.id}\t{codon_num}\t{codon_pos}\n'
                        fh.write(row)
                        total_sites +=1
            # Handle genes in reverse strand
            else:
                # Temporary data for sorting the output
                rows = dict()
                # In this case, we want to traverse the exons in the opposite order
                for cds in transcript.exons_l[::-1]:
                    # We also navigate the positions in the reverse order
                    for i, pos in enumerate(range(cds.stop, cds.start-1, -1)):
                        codon_pos = ((i-cds.phase)%3)+1
                        if codon_pos == 1:
                            codon_num+=1
                        row = f'{transcript.chrom}\t{pos}\t{transcript.id}\t{cds.id}\t{codon_num}\t{codon_pos}\n'
                        # Don't save this yes, so they can be in the right order in the file
                        # Instead, save in the temporary dictionary
                        rows[pos] = row
                        total_sites +=1
                # Once all the negative strand positions have been parsed, then save in the right order
                for pos in sorted(rows.keys()):
                    fh.write(rows[pos])

    print(f'    Calculated codon positions for {total_sites:,} coding sites.', flush=True)

def main():
    print(f'{PROG} started on {date()} {time()}.')
    args = parse_args()
    # Extract the CDS annotations from the GFF
    annotations = load_cds_from_gff(args.gff)
    # Calculate the codon positions
    calculate_codon_positions(annotations, args.out_dir)
    print(f'\n{PROG} finished on {date()} {time()}.')

# Run Code
if __name__ == '__main__':
    main()
