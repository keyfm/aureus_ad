# read subject files with sample IDs (input)
# re-build fa file names 
# merge per subject and bgzip
# run spades on newly generated file (output)

import argparse,sys,gzip,os
from Bio import SeqIO
import statistics
# import subprocess # needed for system calls


''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
    Py Script integrated in de-novo genome assembly snakemake routine.
    Extracts summary statistics for de-novo assembly
                                ''',
                                epilog="Questions or comments? --> fkey@mit.edu")
parser.add_argument("-s", dest='samples', help="File or sample ID string. Input file per subject including sample-IDs validated by kraken. Alternatively, if single samples used for assembly just provide sample ID.")
parser.add_argument("-f", dest='fastq', help="Input FQ of all validated samples per subject. NOTE: bgzip!", type=argparse.FileType('rt'))
parser.add_argument("-c", dest='contig', help="contigs.fa SPAdes output file", type=argparse.FileType('rt'))
parser.add_argument("-a", dest='annotation', help="Summary annotation txt output from prokka", type=argparse.FileType('rt'))
parser.add_argument("-o", "--out", dest="out", help="Output table. If not specified stdout", type=argparse.FileType('w'), default=sys.stdout)
args = parser.parse_args()


''' FUNCTIONS'''
def get_numSample(samples_file):
    ''' count lines samples file '''
    with open(samples_file) as file:
        num_lines = sum(1 for line in file)
    return num_lines

def get_readPairCount(fastq_bgzip):
    ''' read in froward bgzip FQ file...report pair count'''
    with gzip.open(fastq_bgzip,"r") as file:
        num_lines = sum(1 for line in file)
    num_reads = int(num_lines/4)
    return num_reads

def fasta_parse(fasta_file):
    ''' get relevant assembly stats '''
    ls_lgt = []
    ls_cov = []
    for fasta in fasta_file:
        hdr = fasta.id
        hdr = hdr.split('_')
        if int(hdr[3]) > 500:
            # print(hdr[3])
            ls_lgt.append(int(hdr[3]))
            ls_cov.append(float(hdr[5]))
    return [ls_lgt,ls_cov]

def mean_sd_cov(ls_len,ls_cov):
    ''' build vector of 10000 with proportional cov values to estimate cov/SD (avoid building full-genome-length vector)'''
    ls_ratio = []
    for i in ls_len:
        ls_ratio.append( round( i/sum(ls_len)*10000 , 0))
    ls_cov_dis_estimate = []
    # pretty awful nested loop to create a vector of ~10k values representing proportional to contig length the contig-specific cov values
    for i in range(len(ls_ratio)):
        for j in range(int(ls_ratio[i])):
            ls_cov_dis_estimate.append(ls_cov[i])
    # use 10k list representation of cov values to calc SumStats
    mean_cov_estimate = round( sum(ls_cov_dis_estimate) / len(ls_cov_dis_estimate) , 2)
    cov_sd = round(statistics.stdev(ls_cov_dis_estimate),4)
    return [mean_cov_estimate , cov_sd]

def get_n50_l50(ls_len):
    ''' calculate N50/L50 stats '''
    full_len = sum(ls_len)
    half_len = round(full_len/2,0)
    cum_ctr = 0
    l50 = 0
    for i in sorted(ls_len,reverse=True): # reverse just to make sure descending order
        cum_ctr = cum_ctr + i
        l50 = l50 + 1
        if cum_ctr >= half_len:
            n50 = i
            break
    return [n50,l50]

def get_subjectID(filename):
    ''' get subject identifier from file name supplied to py command'''
    subject_identifier = filename.split("/")[-1].split("_")[0].replace('subject','')
    return subject_identifier


def get_annotation_nums(file):
    res_ls = [0,0,0,0] # placeholder for CDS, tRNA, rRNA, tmRNA count
    for line in file:
        line = line.strip('\n')
        line = line.split(' ')
        if line[0].startswith('CDS:'):
            res_ls[0] = line[1]
        if line[0].startswith('tRNA:'):
            res_ls[1] = line[1]
        if line[0].startswith('rRNA:'):
            res_ls[2] = line[1]
        if line[0].startswith('tmRNA:'):
            res_ls[3] = line[1]
    return res_ls
        


''' MAIN '''

if __name__ == "__main__":

    # readinput
    samplesFileOrString = args.samples
    FQgz = args.fastq
    contig = args.contig
    contigs = SeqIO.parse(open(contig.name),'fasta')


    # get sumary stats
    if os.path.isfile(samplesFileOrString):
        sampleCount = get_numSample(samplesFileOrString)
        subjectID = get_subjectID(samplesFileOrString)
    else:
        sampleCount = 1 # assume only one single sample provided
        subjectID = samplesFileOrString
    fqReadCount = get_readPairCount(FQgz.name)
    list_length,list_cov = fasta_parse(contigs) # for contigs > 500bp: extract length and cov as separate list
    assLength = sum(list_length)
    numcontig = len(list_length)
    maxContig = max(list_length)
    meanContigLen = int(round( sum(list_length) / len(list_length) , 0))
    medianContigLen = statistics.median(list_length) 
    meanCov,covSD =  mean_sd_cov(list_length,list_cov)
    n50,l50 = get_n50_l50(list_length)
    [cds,tRNA,rRNA,tmRNA] = get_annotation_nums(args.annotation)
    
    # write output
    args.out.write('\t'.join(['SampleID','NumSamplesValid','NumReadsSpadesInput','AssemblyLength','ContigCount','MaxContig','N50','L50','MeanContigLength','MedianContigLength','MeanCov','SDcov','CDS','tRNA','rRNA','tmRNA']) + '\n')
    args.out.write('\t'.join(map(str,[subjectID,sampleCount,fqReadCount,assLength,numcontig,maxContig,n50,l50,meanContigLen,medianContigLen,meanCov,covSD,cds,tRNA,rRNA,tmRNA])) + '\n')
    args.out.close()

exit()
