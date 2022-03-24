# read subject files with sample IDs (input)
# re-build fa file names 
# merge per subject and bgzip
# run spades on newly generated file (output)

import argparse,sys
import subprocess # needed for system calls


''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
    Py Script integrated in de-novo genome assembly snakemake routine.
    Step1: Merge previously created and kraken-validated fq.gz files (each contains 250k)
    Step2: Run SPAdes in careful mode.
    Note: Input file has to follow naming scheme in order to allow the script to extract the subject identifier (see def run_spades : [folder]/[samplesPerSubject]/samples[subjectID]_[string.txt])
                                ''',
                                epilog="Questions or comments? --> fkey@mit.edu")
parser.add_argument("-i", dest='input', help="Input file per subject including sample-IDs validated by kraken", type=argparse.FileType('rt'))
parser.add_argument('-t', dest='threads',help="Number of threads",type=int,default=1)
parser.add_argument('-s', dest='subjectid',help="Subject ID (Snakemake wildcard!)",type=str)
parser.add_argument('-e', dest='exeSpades',help="Path to spades.py",type=str)
args = parser.parse_args()


''' FUNCTIONS'''

def build_sample_file_list(file):
    ''' read subject file with sample IDs and transform to files (based on known structure) '''
    # with open(file, "r") as ins:
    ls_fwd = []
    ls_rev = []
    for line in file:
        line = line.rstrip('\n')
        ls_fwd.append("../kraken2/0-tmp/"+line+"_1.fq.gz")
        ls_rev.append("../kraken2/0-tmp/"+line+"_2.fq.gz")
    return [ls_fwd,ls_rev]

def merge_fq(ls_sample_file,subject):
    ''' take all sample-files and merge them using system cat call (rev & fwd) '''
    input_file1_string = ' '.join(ls_sample_file[0])
    input_file2_string = ' '.join(ls_sample_file[1])
    outfile1 = "../kraken2/0-tmp/in1_spades_"+subject+".fq.gz"
    outfile2 = "../kraken2/0-tmp/in2_spades_"+subject+".fq.gz"
    subprocess.run("zcat " + input_file1_string + " |bgzip > " + outfile1, shell=True)
    subprocess.run("zcat " + input_file2_string + " |bgzip > " + outfile2, shell=True)
    return [outfile1,outfile2]


def run_spades(executable_spades,ls_merge_fq, subject, threads):
    ''' run spades. output to 3-spades/subject_[]/ '''
    outfolder = "3-spades/subject_"+subject+"/"
    subprocess.run("mkdir -p " + outfolder , shell=True)
    subprocess.run(executable_spades + " --careful -t " + str(threads) + " -1 " + ls_merge_fq[0] + " -2 " + ls_merge_fq[1] + " -o " + outfolder , shell=True)


''' MAIN '''
if __name__ == "__main__":

    # infile = '3-spades/samplesPerSubject/subject10_samples.txt' # from cammand line
    infile = args.input
    file_names = build_sample_file_list(infile)
    subjectID = args.subjectid 
    threads = args.threads
    executable_spades = args.exeSpades
    outfileLs = merge_fq(file_names,subjectID)
    run_spades(executable_spades,outfileLs,subjectID,threads)
    sys.exit()


