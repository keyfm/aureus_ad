# read multifasta generated with spades and output contigs longer than xb

import argparse,sys
from Bio import SeqIO


''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
    Py Script to read multifasta generated with spades and output contigs longer than xb
                                ''',
                                epilog="Questions or comments? --> fkey@mit.edu")
parser.add_argument("-f", dest='fasta', help="multifasta (contigs) SPAdes output file", type=argparse.FileType('rt'))
parser.add_argument('-l', dest='length',help="Min length per contig",type=int,default=1)
parser.add_argument("-o", "--out", dest="out", help="Output fasta. If not specified stdout", type=argparse.FileType('w'), default=sys.stdout)
args = parser.parse_args()


def fasta_parse(fasta_file,min_length):
    ''' get relevant assembly stats '''
    rec_list = []
    for fasta in fasta_file:
        hdr = fasta.id
        hdr = hdr.split('_')
        if int(hdr[3]) >= min_length:
            rec_list.append(fasta)
    return rec_list


if __name__ == "__main__":
    contig = args.fasta
    contigs = SeqIO.parse(open(contig.name),'fasta')
    seqIOrec_passed = fasta_parse(contigs , args.length)
    SeqIO.write(seqIOrec_passed, args.out, "fasta")
