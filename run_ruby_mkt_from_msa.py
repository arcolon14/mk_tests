#!/usr/bin/env python3
import sys, os, argparse, subprocess
from statistics import mean
from statistics import stdev
from datetime import datetime
from warnings import warn
PROG = sys.argv[0].split('/')[-1]
DESC = """Run AK's custom Ruby script for MK-Test anslysis."""

def parse_args():
    '''Set and verify command line options.'''
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-l', '--sco-list', required=True, 
                   help='(str) Path to the file containing list of single-copy orthologs.')
    p.add_argument('-m', '--msa-dir', required=True,
                   help='(str) Path to directory containing MSA files.')
    p.add_argument('-g', '--outgroup-id', required=True,
                   help='(str) ID of the outgroup species in the MSA.')
    p.add_argument('-e', '--exe-dir', required=False, default='.',
                   help='(str) Path to the Ruby `mkTest.rb` executable [default=.].')
    p.add_argument('-o', '--out-dir', required=False, default='.',
                   help='(str) Path to output directory [default=.].')
    # Check inputs
    args = p.parse_args()
    assert os.path.exists(args.sco_list)
    assert os.path.exists(args.out_dir)
    args.out_dir = args.out_dir.rstrip('/')
    assert os.path.exists(args.msa_dir)
    assert os.path.exists(args.exe_dir)
    assert os.path.exists(f'{args.exe_dir}/mkTest.rb')
    return args

def date() -> str:
    '''Print the current date in YYYY-MM-DD format.'''
    return datetime.now().strftime("%Y-%m-%d")

def time() -> str:
    '''Print the current time in HH:MM:SS format.'''
    return datetime.now().strftime("%H:%M:%S")

def load_sco_list(sco_list_f:str)->list:
    '''
    Load the file containing the ID of single-copy orthologs.
    Args:
        sco_list_f: (str) Path to file containing SCO IDs.
    Returns:
        scos: (list) List of SCO IDs.
    '''
    scos = list()
    print('\nParsing list of single-copy orthologs...', flush=True)
    with open(sco_list_f) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            line = line.strip('\n')
            if len(line) == 0:
                continue
            if not line.startswith('N0.HOG'):
                sys.exit('Error: Single-copy ortholog file must contain list of orthofinder N0 orthogroups.')
            scos.append(line)
    # Make checks for duplicate elements
    scos = set(scos)
    scos = list(scos)
    # Report and Return
    print(f'    Extracted {len(scos):,} single-copy orthologs from input file.', flush=True)
    return scos

class mkTestResults:
    '''
    Store the results for an mkTest run.
    '''
    def __init__(self,
                 ortho_id:str,
                 gene_id:str='NA',
                 aln_len:int=-1,
                 aaFix:int=-1,
                 aaPoly:int=-1,
                 silFix:int=-1,
                 silPoly:int=-1,
                 p_val:float=-1.0,
                 len_dist=list(),
                 frameshifts:bool=False,
                 passed:bool=True,
                 failure_reason:str='passed'):
        self.orthoID = ortho_id
        self.geneID  = gene_id
        self.alnLen  = aln_len
        self.aaFix   = aaFix
        self.aaPoly  = aaPoly
        self.silFix  = silFix
        self.silPoly = silPoly
        self.pVal    = p_val
        self.lenDist = len_dist
        self.fshifts = frameshifts
        self.passed  = passed
        self.whyfail = failure_reason
    def __str__(self):
        if len(self.lenDist) > 0:
            mean_len = mean(self.lenDist)
            sd_len = stdev(self.lenDist)
        else:
            mean_len = -1
            sd_len = -1
        return f'{self.orthoID}\t{self.geneID}\t{self.alnLen}\t{mean_len}\t{sd_len}\t{self.aaFix}\t{self.aaPoly}\t{self.silFix}\t{self.silPoly}\t{self.pVal}\t{self.whyfail}'
    def write_row(self):
        if len(self.lenDist) > 0:
            mean_len = mean(self.lenDist)
            sd_len = stdev(self.lenDist)
        else:
            mean_len = -1
            sd_len = -1
        return f'{self.orthoID}\t{self.geneID}\t{self.alnLen}\t{mean_len:0.3f}\t{sd_len:0.3f}\t{self.aaFix}\t{self.aaPoly}\t{self.silFix}\t{self.silPoly}\t{self.pVal:0.8g}\t{self.whyfail}\n'

def read_msa(msa_fa_f:str)->dict:
    '''
    Parse an MSA FASTA and extract the aligned sequences.
    Args:
        msa_fa_f: (str) Path to MSA FASTA file.
    Returns:
        msa_seqs: (dict) Dictionary with the aligned sequences.
    '''
    msa_seqs = dict()
    with open(msa_fa_f) as fh:
        name = None
        seq = ''
        for line in fh:
            line = line.strip('\n')
            if len(line) == 0 or line.startswith('#'):
                continue
            if line.startswith('>'):
                # When reading the next sequence in the FASTA
                if name is not None:
                    msa_seqs[name] = seq
                # Store the ID and prepare for the next sequence
                name = line.lstrip('>')
                seq = ''
            else:
                # If it is a sequence...
                seq += line.upper()
        # At the end of the file, process the remaining sequence.
        msa_seqs[name] = seq
    return msa_seqs

def run_mktest_cmd(result:mkTestResults, ingroup_msa:str, outgroup_msa:str, exe_dir:str)->list:
    '''
    Run the `ruby mkTest.rb` command using the target ingroup and 
    outgroup alignments.
    Args:
        result: (mkTestResults) Object to store mkTest results.
        ingroup_msa: (str) Path to MSA with the ingroup sequences.
        outgroup_msa: (str) Path to MSA with the outgroup sequences.
        exe_dir: (str) Path to directory contain the ruby executable.
    Returns:
        outputs: (list) Outputs of the command.
    '''
    assert isinstance(result, mkTestResults)
    # First, move to the directory of the ruby script
    os.chdir(exe_dir)
    cmd = ['ruby', 'mkTest.rb',
           ingroup_msa,
           outgroup_msa]
    cmd_str = ' '.join(cmd)
    # Then, run the command
    process = None
    # try:
    #     process = subprocess.run(cmd, capture_output=True, check=True, text=True)
    # except subprocess.CalledProcessError as e:
    #     warn(f'Warning: the command\n    {cmd_str}\nfinished with a non-zero status ({e}).\n')
    process = subprocess.run(cmd, capture_output=True, text=True)
    e = process.returncode
    if e != 0:
        msg = f'Warning: the command\n    {cmd_str}\nfinished with a non-zero status ({e}).\n'
        warn(msg)
        result.passed = False
        result.whyfail = 'cmdFailed'
        return result
    stdout = process.stdout
    for line in stdout.strip('\n').split('\n'):
        # Output should look like:
        # aaFix  aaPoly  silFix  silPoly  FET_p_val
        # 6	     1       8       8        1.76e-01
        if len(line) == 0:
            continue
        if line.startswith('aaFix'):
            continue
        fields = line.split('\t')
        # Store as the results class
        result.aaFix   = int(fields[0])
        result.aaPoly  = int(fields[1])
        result.silFix  = int(fields[2])
        result.silPoly = int(fields[3])
        result.pVal    = float(fields[4])
    return result

def split_msa(sco:str, aligned_sequences:dict, outgroup_id:str, out_dir:str, fa_line_width:int=60)->tuple:
    '''
    Split the input MSA into ingroup/outgroup sequences.
    Args:
        sco: (str) Single copy ortholog ID.
        aligned_sequences: (dict) Dictionary with the aligned sequences.
        outgroup_id: (str) ID of ourgroup taxon in the MSA.
        out_dir: (str) Output directory.
        fa_line_width: (int) Character width in output FASTA [default=60]
    Return:
        tuple:
            ingroup_msa: (str) Path to MSA with the ingroup sequences.
            outgroup_msa: (str) Path to MSA with the outgroup sequences.
    '''
    # Prepare the two outputs
    ingroup_msa  = f'{out_dir}/{sco}_ingroup.msa.fa'
    ingroup_fh    = open(ingroup_msa, 'w')
    outgroup_msa = f'{out_dir}/{sco}_outgroup.msa.fa'
    outgroup_fh   = open(outgroup_msa, 'w')
    # Parse the MSA sequences, identify the outgroup, and
    # save to the corresponding file.
    # Some checks:
    outgroup_seen = 0
    ingroup_seen = 0
    for seq_id in aligned_sequences:
        sequence = aligned_sequences[seq_id]
        out_file = ingroup_fh
        # See if outgroup ID is in the FASTA header
        if outgroup_id in seq_id:
            outgroup_seen += 1
            out_file = outgroup_fh
        else:
            ingroup_seen += 1
        # Once the target file has been selected, save sequence.
        out_file.write(f'>{seq_id}\n')
        # Wrap the sequence lines up to `fa_line_width` characters
        for start in range(0, len(sequence), fa_line_width):
            seq_line = sequence[start:(start+fa_line_width)]
            out_file.write(f'{seq_line}\n')
    # Make some checks:
    if outgroup_seen < 1:
        sys.exit(f'Error: {outgroup_id} outgroup not found for ortholog {sco}.')
    if ingroup_seen < 1:
        sys.exit(f'Error: ingroup sequences (non-{outgroup_id}) not found for ortholog {sco}.')
    if outgroup_seen > ingroup_seen:
        warn(f'Warning: More outgroup ({outgroup_seen}) than ingroup ({ingroup_seen}) sequences seen for ortholog {sco}.')
    ingroup_fh.close()
    outgroup_fh.close()
    return ingroup_msa, outgroup_msa

def process_ortholog(sco:str, outgroup_id:str, 
                     msa_dir:str, exe_dir:str,
                     out_dir:str)->mkTestResults:
    '''
    Process a target single-copy ortholog: Split alignment, run MKtest, 
    and parse output.
    Args:
        sco: (str) Single copy ortholog ID.
        outgroup_id: (str) ID of ourgroup taxon in the MSA.
        msa_dir: (str) Path to directory containing MSAs.
        exe_dir: (str) Path to directory contain the ruby executable.
        out_dir: (str) Output directory.
    Returns:
        result: (mkTestResults) Object of the mkTest script results.
    '''
    result = mkTestResults(sco)
    # 1. Read the corresponding MSA
    msa_fa = f'{msa_dir}/{sco}_NT.fa'
    # Some orthologs have no valid MSA, skip them.
    # TODO: Why? MACSE was unable to generate one, it seems.
    if not os.path.exists(msa_fa):
        result.passed = False
        result.whyfail = 'noMSA'
        return result
    aligned_seq = read_msa(msa_fa)
    # 2. Inspect the sequence:
    #    Get the target gene and alignment length.
    #    Look for frameshifts and other issues
    # For the lengths of the aligned seqs
    result.lenDist = [ 0 for _ in range(len(aligned_seq)) ]
    for i, seq_id in enumerate(aligned_seq):
        sequence = aligned_seq[seq_id]
        for s in sequence:
            # Tally the non-gap length
            if s.upper() in 'ACGTN':
                result.lenDist[i] += 1
            # Check for frameshifts
            elif s.upper() == '!':
                result.fshifts = True
        # If not the outgroup, extract length and 
        # gene information
        if outgroup_id not in aligned_seq:
            result.alnLen = len(sequence)
            # Note: this works for our default input, but 
            # will likely fail for other cases.
            # TODO: Make more robust?
            gene_ID = seq_id.split('_')[0]
            if not gene_ID.startswith('mrna'):
                warn(f'Warning: Non-default FASTA headers found for ortholog {sco}.')
            result.geneID = gene_ID
    # Skip sequences with frameshifts
    if result.fshifts:
        result.passed = False
        result.whyfail = 'frameshifts'
        return result
    # 3. Split the MSAs
    ingroup_msa, outgroup_msa = split_msa(sco, aligned_seq, outgroup_id, out_dir)
    # 4. Run the mkTest.rb command.
    result = run_mktest_cmd(result, ingroup_msa, outgroup_msa, exe_dir)
    return result

def process_all_orthologs(scos:list, outgroup_id:str, 
                          msa_dir:str, exe_dir:str, 
                          out_dir:str)->None:
    '''
    Process all single-copy orthologs: Split alignment, run MKtest, 
    and parse output.
    Args:
        scos: (list) list of single copy ortholog IDs.
        outgroup_id: (str) ID of ourgroup taxon in the MSA.
        msa_dir: (str) Path to directory containing MSAs.
        exe_dir: (str) Path to directory contain the ruby executable.
        out_dir: (str) Output directory.
    Returns:
        None
    '''
    print('\nProcessing all the single-copy orthologs (splitting alignments, running MK-test, and processing output)...', flush=True)
    n_processed = 0
    # Prepare output directory for the temporary MSAs
    split_dir = f'{out_dir}/split_msa'
    if not os.path.exists(split_dir):
        os.mkdir(split_dir)
    # Process all the data
    with open(f'{out_dir}/rb_mkTest_out.tsv', 'w') as fh:
        header = ['orthologID', 'focalGene', 'alnLen',
                  'avgAlnLenNoGaps', 'sdAlnLenNoGaps',
                  'aaFix', 'aaPoly', 'silFix', 'silPoly',
                  'fetPval', 'notes']
        header = '\t'.join(header)
        fh.write(f'{header}\n')
        # Process the individual orthologs
        for i, sco in enumerate(sorted(scos)):
            # Report a tally of progress
            if i%1000==0 and i>0:
                print(f'    Processing the {i:,}th ortholog.', flush=True)
            # Process a single ortholog
            result = process_ortholog(sco, outgroup_id, 
                                       msa_dir, exe_dir,
                                       split_dir)
            if result.whyfail == 'passed':
                n_processed+=1
            fh.write(result.write_row())
    print(f'    \nSuccesfully processed results for {n_processed:,} single-copy orthologs.', flush=True)

def main():
    print(f'{PROG} started on {date()} {time()}.')
    args = parse_args()
    # Load the single copy orthologs
    scos = load_sco_list(args.sco_list)
    # Process all the data
    process_all_orthologs(scos, args.outgroup_id, args.msa_dir,
                          args.exe_dir, args.out_dir)
    print(f'\n{PROG} finished on {date()} {time()}.')

# Run Code
if __name__ == '__main__':
    main()
