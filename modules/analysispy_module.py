# module accpanying analysis.py
# Module supported and maintaned for analysis of candidate_mutation_table.pickle.gz (created by build_candidate_mutation_table.py; previous matlab-based table not recommended to use due 1-based/0-based conflict!)



import h5py,gzip, os, re, scipy
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from BCBio.GFF import GFFExaminer # pip install bcbio-gff
from BCBio import GFF # might be redundant
import sys
import glob
import pandas as pd
import pickle
from collections import OrderedDict
from Bio import Phylo
from Bio import AlignIO
import subprocess
from pylab import *
import random
from scipy import stats
#from datetime import datetime, date, time
import time,datetime
from matplotlib import rc
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.font_manager import FontProperties
from Bio.Data import CodonTable
import collections
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor #NJ tree
import warnings
import seaborn as sns
from functools import partial # for gff_parse()

def dnaparse2lol(file):
    # Each entry will be element in list. Each element in list is a list of length 2, with header [0] and sequence [1]
    # First header line is skipped
    ctr=0
    fa_lol = [] 
    for line in file:
        ctr += 1
        if ctr == 1:
            continue # skip header
        if ctr > 1:
            line = line.strip().split() # read in reference seq
            fa_lol.append( [line[0],line[1]] )
    return fa_lol

def read_samplesCSV(spls):
    # reads in samples.csv file, format: Path,Sample,ReferenceGenome,ProviderName,Subject
    hdr_check = ['Path', 'Sample', 'ReferenceGenome', 'ProviderName', 'Subject']
    switch = "on"
    file = open(spls, 'r')
    list_path = []
    list_splID = []
    list_providerNames = []
    list_refG = []
    list_patient = []
    for line in file:
        line = line.strip('\n').split(',')
        # Test Header. Note: Even when header wrong code continues (w/ warning), but first line not read.
        if switch == "on":
            if (line == hdr_check):
                print("Passed CSV header check")
            else:
                Warning("CSV did NOT pass header check! Code continues, but first line ignored")
            switch = "off"
            continue
        # build lists
        list_path.append(line[0])
        list_splID.append(line[1])
        list_refG.append(line[2])
        list_providerNames.append(line[3])
        list_patient.append(line[4])
    return [list_path,list_splID,list_refG,list_providerNames,list_patient] # set(list_patient) provides only unique subject IDs

def read_candidate_mutation_table(matlab_file):
    # reads the candidate_mutation_table.m and outputs each dataset as array. 
    # Note: SampleNames == list[]
    f = h5py.File(matlab_file)
    # read data2dict
    arrays = {}
    for k, v in f.items():
        arrays[k] = np.array(v)
    # assign variaus variables. NOTE: import transposes arrays (only x/y but not z axis)
    quals = arrays['Quals'] # transposed. all values negative > keep > TDL changed that later in analysis
    quals = quals.transpose()
    p = arrays['p'] # Note: p is 1-based
    p = p.flatten().astype(np.int64) 
    counts = arrays['counts'] #[level,row,col] > in matlab: level=sample,row=ATGCatgc, col=pos; NOTE: Arolyns PDF 3d object wrong axis!
    if len(counts.shape) != 3: # only prints correct output when used after indel_counter implemented
        print('Attention: counts 3D matrix is not 3D. Skip subject!')
        return [[],[],[],[],[],[],True]
    counts = counts.transpose(0, 2, 1) # transpose: (level==; row>col; col>row)
    in_outgroup = arrays['in_outgroup'].transpose().astype(int) # data in one row for each sample. integer!    
    # SampleNames saved in specific object pointer format which requires this technical loop below to be resolved
    sampleNames = []
    mygroup = f['SampleNames']
    for s in mygroup:
        obj=f[s[0]]
        str1 = ''.join(chr(i) for i in obj[:])
        sampleNames.append(str1)
    if 'indel_counter' in arrays.keys():
        indel_counter = arrays['indel_counter']
        indel_counter = indel_counter.transpose(0, 2, 1) # transpose: (level==; row>col; col>row)
        return [ quals, p, counts, in_outgroup, sampleNames,indel_counter,False]
    else:
        return [ quals, p, counts, in_outgroup, sampleNames]

def read_candidate_mutation_table_pickle_gzip(file_cmt_pickle_gz):
    # read the candidate_mutation_table.pickle.gz files
    with gzip.open(file_cmt_pickle_gz, 'rb') as f:
        cmt = pickle.load(f)
        sampleNames = np.array(cmt[0])
        p = np.array(cmt[1]) # p 0-based index of genomic position (means actual position is p+1)
        counts = np.array(cmt[2])
        quals = np.array(cmt[3])
        in_outgroup = np.array(cmt[4])
        indel_counter = np.array(cmt[5])
    return [quals,p,counts,in_outgroup,sampleNames,indel_counter]


        
def genomestats(REFGENOMEFOLDER):
    # parse ref genome to extract relevant stats
    # accepts genome.fasta or genome.fasta.gz (gzip) in refgenomefolder
    fasta_file = glob.glob(REFGENOMEFOLDER + '/genome.fasta')
    if len(fasta_file) != 1:
        fasta_file_gz = glob.glob(REFGENOMEFOLDER + '/genome.fasta.gz')
        if len(fasta_file_gz) != 1:
            raise ValueError('Either no genome.fasta(.gz) or more than 1 genome.fasta(.gz) file found in ' + REFGENOMEFOLDER)
        else: # genome.fasta.gz
            refgenome = SeqIO.parse(gzip.open(fasta_file_gz[0], "rt"),'fasta')
    else: # genome.fasta
        refgenome = SeqIO.parse(fasta_file[0],'fasta')
    Genomelength = 0
    ChrStarts = []
    ScafNames = []
    for record in refgenome:
        ChrStarts.append(Genomelength) # chr1 starts at 0 in analysis.m
        Genomelength = Genomelength + len(record)
        ScafNames.append(record.id)
    # close file
    #refgenome.close() # biopy update SeqIO has no close attribute anymore.
    # turn to np.arrys!
    ChrStarts = np.asarray(ChrStarts,dtype=int)
    Genomelength = np.asarray(Genomelength,dtype=int)
    ScafNames = np.asarray(ScafNames,dtype=object)
    return [ChrStarts,Genomelength,ScafNames]

def p2chrpos(p, ChrStarts):
    '''# return 2col array with chr and pos on chr
    #p...continous, ignores chr
    #pos: like p, 0-based'''
        
    # get chr and pos-on-chr
    chr = np.ones(len(p),dtype=int)
    if len(ChrStarts) > 1:
        for i in ChrStarts[1:]:
            chr = chr + (p > i) # when (p > i) evaluates 'true' lead to plus 1 in summation. > bcs ChrStarts start with 0...genomestats()
        positions = p - ChrStarts[chr-1] # [chr-1] -1 due to 0based index
        pos = np.column_stack((chr,positions))
    else:
        pos = np.column_stack((chr,p))
    return pos

def extract_outgroup_mutation_positions(REFGENOMEFOLDER,position,CMTpy=True):
    # extracts the ref nucleotide for every position. positions needs to be sorted by chr
    # CMTpy=True: if old matlab build_candidate_mutation.mat used, put flag False. p 1-based correction
    # reads genome.fasta or genome.fasta.gz    
    fasta_file = glob.glob(REFGENOMEFOLDER + '/genome.fasta')
    if len(fasta_file) != 1:
        fasta_file_gz = glob.glob(REFGENOMEFOLDER + '/genome.fasta.gz')
        if len(fasta_file_gz) != 1:
            raise ValueError('Either no genome.fasta(.gz) or more than 1 genome.fasta(.gz) file found in ' + REFGENOMEFOLDER)
        else: # genome.fasta.gz
            refgenome = SeqIO.parse(gzip.open(fasta_file_gz[0], "rt"),'fasta')
    else: # genome.fasta
        refgenome = SeqIO.parse(fasta_file[0],'fasta')
    refnt = np.zeros(position.shape[0],dtype=object)
    pos_counter = 0
    chr_counter = 1
    for record in refgenome:
        # print(record)
        poschr = position[ position[:,0]==chr_counter , 1]
        for sglpos in poschr:
            # print(sglpos)
            if CMTpy:
                refnt[pos_counter] = record.seq[sglpos]._data # when build_candidate_mutation_table.py is used!
            else:
                refnt[pos_counter] = record.seq[sglpos-1]._data # NOTE: uposition based on p, which is 1-based coordinate
            pos_counter += 1
        chr_counter += 1
    return refnt            

def extract_seq_from_fasta(refgenomefolder,chrom,start,end,strand):
    # returns seq object between start and end
    # chrom is 0 based integer! start end 1-based (I think!)
    # NOTE: coordinates relative on chromosome, if more than 1 present in fasta
    # reads genome.fasta or genome.fasta.gz    
    fasta_file = glob.glob(refgenomefolder + '/genome.fasta')
    if len(fasta_file) != 1:
        fasta_file_gz = glob.glob(refgenomefolder + '/genome.fasta.gz')
        if len(fasta_file_gz) != 1:
            raise ValueError('Either no genome.fasta(.gz) or more than 1 genome.fasta(.gz) file found in ' + REFGENOMEFOLDER)
        else: # genome.fasta.gz
            refgenome = SeqIO.parse(gzip.open(fasta_file_gz[0], "rt"),'fasta')
    else: # genome.fasta
        refgenome = SeqIO.parse(fasta_file[0],'fasta')
    for c,record in enumerate(refgenome):
        if c == chrom:
            if strand == 1:
                return record.seq[int(start)-1:int(end)]
            elif strand == -1:
                return record.seq[int(start)-1:int(end)].reverse_complement()
            else:
                print('Unresolved strandedness.')


def div_major_allele_freq(cnts):
    # define matrices with major allele freq; major Nucl. (maNT); minor Nucl.; minor AF
    
    c=cnts[:,0:4,:]+cnts[:,4:8,:]; # flatten frw and rev ATCG counts    

    sorted_arr = np.sort(c,axis=1) # return sorted_arr matrix; axis 0 in 3d array == sort by col
    sortedpositions = np.argsort(c,axis=1) # return matrix indices of sort;axis 0 in 3d array == sort by col
    
    maxcount = sorted_arr[:,3:4:,:] # get allele counts for major allele (4th row); weird "3:4:" indexing required to maintain 3d structure
    minorcount = sorted_arr[:,2:3:,:] # get allele counts for first minor allele (3rd row); tri/quadro-allelic ignored; weird "2:3:" indexing required to maintain 3d structure and get 3rd row
    
    with np.errstate(divide='ignore', invalid='ignore'): # suppress warning for division by 0
        maf = maxcount / sorted_arr.sum(axis=1,keepdims=True)
        minorAF = minorcount / sorted_arr.sum(axis=1,keepdims=True)

    maf = np.squeeze(maf,axis=1) # turn 2d; axis=1 to keep 2d structure when only one position!
    maf[ np.isnan(maf) ] = 0 # set to 0 to indicate no data
    minorAF = np.squeeze(minorAF,axis=1) 
    minorAF[ np.isnan(minorAF) ] = 0 # set to 0 to indicate no data/no minor AF
    
    majorNT = np.squeeze(sortedpositions[:,3:4:,:],axis=1) # index position in sortedpositions represents allele position ATCG; axis=1 to keep 2d structure when only one position!
    minorNT = np.squeeze(sortedpositions[:,2:3:,:],axis=1) # -"-

    # Note: If counts for all bases are zero, then sort won't change the order
    # (since there is nothing to sort), thus majorNT.minorNT put to 4 (ie NA) using maf (REMEMBER: minorAF==0 is a value!)

    majorNT[maf==0] = 4
    minorNT[maf==0] = 4    
    
    return [maf.transpose(), majorNT.transpose(), minorNT.transpose(), minorAF.transpose()]


def major_allele(arr):
    ''' returns 1-dimensional array of size arr.shape[0] with the major allele index [0:3] or NA [4] (if ambigous)'''
    # NOTE: if major NT ambigous (multiple alleles with same number of occurence I report allele with lower numeric value. could be improved as it might lead to negligible (imo) bias). Also loop could prob be removed.
    # input 2d arr (2nd dim can be 1!) with nucleotide indices (0:4)
    nonNA_out = (arr != 4)
    out_idx = []
    for i,row in enumerate(arr):
        if np.any(nonNA_out[i,:]): # any majorNT found
            row = row[nonNA_out[i,:]]
            row_ct = np.unique(row,return_counts=True)
            idx = np.where(row_ct[1] == np.max(row_ct[1]) )[0]
            out_idx.append(row_ct[0][idx][0]) # if multiple outgroup alleles same count I take the first. Could be changed to NA or add more outgroup samples for refined inference.
        else: # no majorNT
            out_idx.append(4)
    out_idx = np.array(out_idx)
    print(np.sum(out_idx == 4),'/', out_idx.size ,'elements of p have no major allele (ie. 4)!'  )
    return out_idx



def ana_mutation_quality(Calls,Quals):
    # This functions aims at providing the FQ value for every SNP position 
    # Across all pairwise different allele calls, it reports the best FQ value among the minimum FQ values per pair
    # NOTE: This function requires some more efficient coding!
    [Nmuts, NStrain] = Calls.shape ;
    MutQual = np.zeros((Nmuts,1)) ; 
    MutQualIsolates = np.zeros((Nmuts,2)); 
    
    # generate template index array to sort out strains gave rise to reported FQ values
    s_template=np.zeros( (len(Calls[0,:]),len(Calls[0,:])) ,dtype=object)
    for i in range(s_template.shape[0]):
        for j in range(s_template.shape[1]):
            s_template[i,j] = str(i)+"_"+str(j)

    for k in range(Nmuts):
        if len(np.unique(np.append(Calls[k,:], 4))) <= 2: # if there is only one type of non-N (4) call, skip this location
            MutQual[k] = np.nan ;
            MutQualIsolates[k,:] = 0; 
        else:
            c = Calls[k,:] ; c1 = np.tile(c,(c.shape[0],1)); c2 = c1.transpose() # extract all alleles for pos k and build 2d matrix and a transposed version to make pairwise comparison
            q = Quals[k,:] ; q1 = np.tile(q,(q.shape[0],1)); q2 = q1.transpose() # -"-
            g = np.all((c1 != c2 , c1 != 4 , c2 != 4) ,axis=0 )  # no data ==4; boolean matrix identifying find pairs of samples where calls disagree (and are not N) at this position
            #positive_pos = find(g); # numpy has no find; only numpy where, which does not flatten 2d array that way
            # get MutQual + logical index for where this occurred
            MutQual[k] = np.max(np.minimum(q1[g],q2[g])) # np.max(np.minimum(q1[g],q2[g])) gives lower qual for each disagreeing pair of calls, we then find the best of these; NOTE: np.max > max value in array; np.maximum max element when comparing two arryas
            MutQualIndex = np.argmax(np.minimum(q1[g],q2[g])) # return index of first encountered maximum!
            # get strain ID of reorted pair (sample number)
            s = s_template
            strainPairIdx = s[g][MutQualIndex]
            MutQualIsolates[k,:] = [strainPairIdx.split("_")[0], strainPairIdx.split("_")[1]]
            
    return [MutQual,MutQualIsolates]



def parse_gff(REFGENOMEFOLDER,ScafNames,ortholog_info_series=pd.Series(),forceReDo=False):
    # parse gff file (tested with version3) with genome annotation (gff or gff.gz)
    # requires one gff file in REFGENOMEFOLDER, which is detected automatically
    # provide ScafNames to maintain same order of contigs/chr/scaffolds in output list as in previously generated ref-genome-based variables, eg. contig_positions,ChrStarts, GenomeLength
    # some columns potentially contain multiple entries. Those are separated by ";"
    # no data is always reported as "."
    # https://biopython.org/wiki/GFF_Parsing 
    # NOTE: annotation is read from gff file!
    # only execute if dataframes not yet made
    if os.path.isfile(REFGENOMEFOLDER+"/annotation_genes.pandas.py.pk1") and (forceReDo == False):
        afile = open(REFGENOMEFOLDER+"/annotation_genes.pandas.py.pk1", 'rb')
        list_of_dataframes = pickle.load(afile)
        afile.close()
        return list_of_dataframes
    else: # search for gff or gff.gz
        examiner = GFF.GFFExaminer()        
        gff_file = glob.glob(REFGENOMEFOLDER + '/*.gff*') # search for gff or gff.gz or gff3
        if len(gff_file) != 1:
            raise ValueError('Either no gff(.gz) file or more than 1 *gff(.gz) file found in ' + REFGENOMEFOLDER)
        if gff_file[0][-2:] == 'gz':
            encoding = 'gzip'
        else:
            encoding = 'unzip'
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open # define opening procedure gzip.open or open
        with _open(gff_file[0]) as gff_handle:
            possible_limits = examiner.available_limits(gff_handle)
        # start processing
        chromosomes = possible_limits["gff_id"].keys() # list of tuples containing contig ID; do not loop over those, as I only loop over contigs observed (ScafNames)...might need a fix if !=
        tagnumber_counter = 0 # unique, continous numerical id for all features across all chromosomes
        list_of_dataframes = [] # output: each element contains a dtaframe w/ all annotation info for chr. ordered as in ScafNames.
        # define gff_type to extracted. Use all but gene and region. gene annotated in NCBI but lacks most info that extra CDS has. Region just describes chromosome extensions
        limits = dict(gff_type = [i[0] for i in possible_limits['gff_type'].keys() if i[0] != 'gene' and i[0] != 'region'] )
        # loop over chr with established order
        for chrom in ScafNames:
            limits["gff_id"] = [chrom]
            # limits = dict(gff_id=[chrom])
            with _open(gff_file[0]) as gff_handle:
                for rec in GFF.parse(gff_handle, limit_info=limits):
                    # for loop over every chr but only defined [limits] has data. Memory efficient!
                    if rec.id == chrom:
                        # if chr has any feature build list of dicts and append to list_of_dataframes, else append empty dataframe
                        if len(rec.features) > 0:
                            # test if seq object part of gff (prokka-based yes, but NCBI-based no >> then load ref genome.fasta)
                            if len(rec.seq) == rec.seq.count('?'):
                                for seq_record in SeqIO.parse(REFGENOMEFOLDER+"/genome.fasta", "fasta"):
                                    if seq_record.id == rec.id:
                                        rec.seq = seq_record.seq
                                if len(rec.seq) == rec.seq.count('?'): # test if succesful
                                    print('Warning: No reference genome found that matches chromosome:' + rec.id)
                            print(rec.id)
                            lod_genes = [] # list-of-dictionary. Easy to turn to pd.dataframe
                            for gene_feature in rec.features:
                                gene_dict = {}
                                tagnumber_counter += 1
                                gene_dict['type'] = gene_feature.type
                                gene_dict['locustag'] = gene_feature.id
                                # add ortholog info if locustag (eg. repeat region has none)
                                #print("gene_feature.id   "+gene_feature.id)
                                if gene_feature.id != "" and gene_feature.type == 'CDS' and not ortholog_info_series.empty:
                                    gene_dict['orthologtag'] = ortholog_info_series[ortholog_info_series.str.findall(gene_feature.id).str.len() == 1].index[0]
                                #print(rec.id+"   "+gene_dict['locustag'])
                                if 'gene' in gene_feature.qualifiers.keys():
                                    gene_dict['gene'] = ";".join(gene_feature.qualifiers['gene'])
                                else:
                                    gene_dict['gene'] = "." # add "." instead of []
                                gene_dict['type'] = gene_feature.type
                                gene_dict['locustag'] = gene_feature.id
                                if gene_dict['type'] == "CDS" or gene_dict['type'] == "gene":
                                    gene_dict['tagnumber'] = tagnumber_counter
                                else:
                                    gene_dict['tagnumber'] = 0
                                if 'product' in gene_feature.qualifiers.keys():
                                    gene_dict['product'] = ";".join(gene_feature.qualifiers['product']) 
                                else:
                                    gene_dict['product'] = "."
                                if 'protein_id' in gene_feature.qualifiers.keys():
                                    gene_dict['protein_id'] = gene_feature.qualifiers['protein_id']
                                else:
                                    gene_dict['protein_id'] = "."
                                if "Dbxref" in gene_feature.qualifiers.keys():
                                    gene_dict['db_xref'] = ";".join(gene_feature.qualifiers['Dbxref'])
                                else:
                                    gene_dict['db_xref'] = "."
                                if "note" in gene_feature.qualifiers.keys():
                                    gene_dict['note'] = ";".join(gene_feature.qualifiers['note'])
                                elif "Note" in gene_feature.qualifiers.keys():
                                    gene_dict['note'] = ";".join(gene_feature.qualifiers['Note'])
                                else:
                                    gene_dict['note'] = "."
                                if 'phase' in gene_feature.qualifiers.keys():
                                    gene_dict['codon_start'] = int(gene_feature.qualifiers['phase'][0])
                                else:
                                    gene_dict['codon_start'] = "."
                                gene_dict['indices'] = [gene_feature.location.start.position,gene_feature.location.end.position]
                                gene_dict['loc1'] = gene_feature.location.start.position # automatically 0-based
                                gene_dict['loc2'] = gene_feature.location.end.position
                                gene_dict['location'] = gene_feature.location
                                gene_dict['strand'] = gene_feature.location.strand
                                dna_seq = rec.seq[gene_feature.location.start:gene_feature.location.end]
                                if gene_dict['strand'] == 1:
                                    gene_dict['sequence'] = dna_seq
                                elif gene_dict['strand'] == -1:
                                    gene_dict['sequence'] = dna_seq.reverse_complement()
                                else:
                                    gene_dict['sequence'] = dna_seq # eg. repeat region
                                # translation, add info where codon starts if info was available. Usually it starts at 0
                                if isinstance( gene_dict['codon_start'] , int):
                                    sequence2translate = gene_dict['sequence'][gene_dict['codon_start']:]
                                    gene_dict['translation'] = sequence2translate.translate(table="Bacterial") # bacterial genetic code GTG is a valid start codon, and while it does normally encode Valine, if used as a start codon it should be translated as methionine. http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:translation
                                elif gene_dict['type'] == "CDS":
                                    sequence2translate = gene_dict['sequence']
                                    gene_dict['translation'] = sequence2translate.translate(table="Bacterial")
                                else:
                                    gene_dict['translation'] = "." # all non-CDS (RNA's or repeat regions) not translated (as those are sometimes also off-frame)
                                lod_genes.append(gene_dict)
                            # after all features parsed > add as dataframe to output list
                            # ATTENTION: SORT pandas DF (annotation not necessarily sorted)
                            df_sort = pd.DataFrame(lod_genes)
                            df_sort = df_sort.sort_values(by=['loc1'])
                            list_of_dataframes.append(df_sort)
                            
                        else:
                            list_of_dataframes.append( pd.DataFrame([{}]) )
        afile = open(REFGENOMEFOLDER+"/annotation_genes.pandas.py.pk1", 'wb')
        pickle.dump(list_of_dataframes, afile)
        afile.close()
        return list_of_dataframes

def idx2nts(calls,missingdata="?"):
    # translate index array to array containing nucleotides
    nucl = np.array(['A','T','C','G',missingdata],dtype=object) # add 5th element --> no data! == index 4
    palette = [0,1,2,3,4] # values present in index-array
    index = np.digitize(calls.ravel(), palette, right=True)
    return nucl[index].reshape(calls.shape)

def nts2idx(nts_array,missingdata="?"):
    # translate nucl array to array containing nucleotide indices
    nucl = np.array(['A','T','C','G',missingdata],dtype=object) # add 5th element --> no data! == index 4
    indices = [0,1,2,3,4] # values present in index-array
    for i in indices:
        nts_array[ nts_array == nucl[i] ] = i
    return nts_array.astype(int)


def genomic_position_all(anno_genes_ls,genomelength,chrstarts):
    # get initial placeholder for all genes (differentiate by chr and genic/nongenic)  
    # UPDATE: gene_num > 0 means CDS. changed 'tagnumber' assignment in parse_gff()
    locus_tag_nums = np.ones(genomelength,dtype=float)*0.5 # genome-wide unique CDS tag ('tagnumber'), and intergenic 0.5
    cds_indices = np.ones(genomelength,dtype=float)*0.5 # per chromosome unique; intergenic #chr.5-1; intragenic #gene
    
    for i,this_chr_df in enumerate(anno_genes_ls):
        if this_chr_df.empty: #Skip contigs w/o CDS annotation.
            continue
        
        gene_num = this_chr_df[['tagnumber']].values.flatten() # tdl called cds_num but actually it was not always cds. e.g tRNA; .values.flatten() #1 turn pd.DF to 2D numpy array; #2 turn 2d to 1D 
        genestarts = chrstarts[i] + this_chr_df[['loc1']].values.flatten() ; # .values.flatten() #1 turn pd.DF to 2D numpy array; #2 turn 2d to 1D
        geneends = chrstarts[i] + this_chr_df[['loc2']].values.flatten() ; # .values.flatten() #1 turn pd.DF to 2D numpy array; #2 turn 2d to 1D
        
        for j in range(len(genestarts)-1): # loop over all but last element
            locus_tag_nums[ (genestarts[j]-1):geneends[j] ] = gene_num[j]; # set gene number. Consecutive numbering across all chromosomes; intergenic regions left 0.5

            cds_indices[ (genestarts[j]-1):geneends[j] ] = j+1; # j zero based,but represents chr counter (1-based)
            cds_indices[ (geneends[j]):(genestarts[j+1]-1) ] = j+1+0.5 # intergenic are #chr+0.5; starts with 0.5 prior chr1 which is set when array created
        # mark positions for last gene of chr
        locus_tag_nums[ (genestarts[-1]-1):geneends[-1] ] = gene_num[-1];
        cds_indices[ (genestarts[-1]-1):geneends[-1] ] = len(gene_num); # j zero based,but represents chr counter (1-based)
        # mark trailing cds_intergenic from last gene till end of chromosome (if last chr >> till end of genome)
        if ((i+1) < len(chrstarts)):
            cds_indices[ geneends[-1]:chrstarts[i+1] ] = len(gene_num)+0.5;
        else:
            cds_indices[ geneends[-1]:genomelength ] = len(gene_num)+0.5;
    return [locus_tag_nums,cds_indices]

def annotate_mutations(annotation_genes , p_gp , refnti_gp , ancnti_gp , calls_gp , counts_gp , hasmutation_gp , mutQual_gp, promotersize, ref_genome_folder):
    ''' produce dataframe with annotation info for each goodposition used for tree
    all positions are 1-based! (pos, nt_pos, loc1/2, distance1/2)
    # function combines TDL: annotate_mutations_gb.m and append_annotations.m '''
    # p_gp: genomic position of goodpos
    # p_gp=p[goodpos2use]
    # refnti_gp=refnti_m[goodpos2use,:]
    # calls_gp=calls[goodpos2use,:]
    # counts_gp=counts[:,:,goodpos2use]
    # hasmutation_gp=hasmutation[goodpos2use,:]
    # mutQual_gp = mutQual[goodpos2use,].flatten() 
   # ref_genome_folder=REFGENOMEFOLDER
    [maf_gp, maNT_gp, minorNT_gp, minorAF_gp] = div_major_allele_freq(counts_gp)
    nts = ['A','T','C','G'] #'atcg';
    rc = ['T','A','G','C']  # 'tagc';
    [chrstarts, genomelength,scafnames]= genomestats(ref_genome_folder);
    [locus_tag,chr_tag] = genomic_position_all(annotation_genes, genomelength, chrstarts); #loc_tag numbers all genes across genome (intergenic:0.5); chr_tag: numbers genes per chr (intergenic: #chr+0.5)
    mut_genenum = chr_tag[p_gp]
    loc_tag_p_gp = locus_tag[p_gp]
    pos_chr = p2chrpos(p_gp,chrstarts) # get two-col table chr,pos   
    lod_mutAnno = []
    for i,pos in enumerate(p_gp):
        #print(i,pos)
        mut_annotations = {}
        mut_annotations['gene_num'] = mut_genenum[i]
        mut_annotations['gene_num_global'] = loc_tag_p_gp[i]
        mut_annotations['chr'] = pos_chr[i][0]
        mut_annotations['pos'] = pos_chr[i][1] + 1 # turn 0-based pos into 1-based pos.
        mut_annotations['quals'] = mutQual_gp[i]
        mut_annotations['nts'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_ref'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_anc'] = "." # default. overwritten below if NT defined.
        mut_annotations['nt_alt'] = "." # default. overwritten below if NT defined.
        p_chr = pos_chr[i][0]
        if mut_genenum[i] == int(mut_genenum[i]): # intragenic
            p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])-1] # both -1 necessary bcs list of df and df 0-based
            mut_annotations['product'] = p_anno.loc['product']
            mut_annotations['gene'] = p_anno.loc['gene']
            mut_annotations['protein_id'] = p_anno.loc['protein_id']
            mut_annotations['strand'] = p_anno.loc['strand']
            mut_annotations['loc1'] = p_anno.loc['loc1'] + 1 # 1-based first position
            mut_annotations['loc2'] = p_anno.loc['loc2'] # last position of gene (inclusive)
            mut_annotations['sequence'] = p_anno.loc['sequence']
            mut_annotations['protein_id'] = p_anno.loc['protein_id']
            mut_annotations['note'] = p_anno.loc['note']
            mut_annotations['locustag'] = p_anno.loc['locustag']
            if 'orthologtag' in p_anno:
                mut_annotations['orthologtag'] = p_anno.loc['orthologtag']
            mut_annotations['translation'] = p_anno.loc['translation']
            if p_anno.loc['strand'] == 1: # get position within gene, consider strandedness
                mut_annotations['nt_pos'] = mut_annotations['pos'] - mut_annotations['loc1']+1 # pos/loc1 1-based. nt_pos 1-based. +1 to get 1-based nt_pos in gene (checked)!
            elif p_anno.loc['strand'] == -1:
                mut_annotations['nt_pos'] = p_anno.loc['loc2'] - mut_annotations['pos'] +1 # pos/loc2 1-based. +1 to get 1-based nt_pos in gene (checked)!
            else:
                mut_annotations['nt_pos'] = "." # can happen. eg. Crispr
            if mut_annotations['nt_pos'] == ".": # Observed with special 'type's like crispr annotated. rare! leads to no AA inference.
                aan = 9999999
            else:
                aan = int( (mut_annotations['nt_pos']-1 )/3 ) + 1 # get #codon that harbours mutation. 1-based.
                ncn = mut_annotations['nt_pos'] - ((aan-1)*3) # get #nucl within triplett. 1-based
                mut_annotations['aa_pos'] = aan
            codons_ls = []
            aa_ls = []
            if len(mut_annotations['sequence']) >= (aan*3) and mut_annotations['translation'] != ".":
                codon = mut_annotations['sequence'][aan*3-2-1:aan*3] # -1 bcs seq-object 0-based; but positional argument aan not 
                codon = [n for n in codon  ] # turn seq.object to list, seq object does not allow reassignment
                for n in range(4): # test all four nucleotide options for resulting AA change
                    if p_anno.loc['strand'] == 1:
                        codon[ncn-1] = nts[n]
                    else:
                        codon[ncn-1] = rc[n]
                    codonSeqO = Seq( "".join(codon), IUPAC.unambiguous_dna) 
                    codons_ls.append(codonSeqO)
                    aa_ls.append(codonSeqO.translate())
            mut_annotations['codons'] = codons_ls
            mut_annotations['AA'] = [aa[0] for aa in aa_ls]
            # append_annotations intragenic
            if len(mut_annotations['AA']) < 4:
                mut_annotations['type'] = 'U'
            mut_annotations['AA_gt'] = '' # Fill in annotations with whether NS or Syn mutation
            if np.unique(refnti_gp[i,:])[0] != 4: # only if ref defined; [0] ok...see next line comment
                mut_annotations['nt_ref'] = nts[np.unique(refnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                mut_annotations['nts'] = mut_annotations['nt_ref']
                if np.unique(ancnti_gp[i,:])[0] != 4:
                    mut_annotations['nt_anc'] = nts[np.unique(ancnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                if len(mut_annotations['AA']) == 4:
                    mut_annotations['AA_gt'] = mut_annotations['AA_gt'] + mut_annotations['AA'][ np.unique(refnti_gp[i,:])[0] ] # the AA list order corresponds to NTS list!
            # extract derived genotype(s) and according AA across all samples
            for j,callidx in enumerate(calls_gp[i,:]):
                #print(str(j)+" "+str(callidx))
                #fixed mutation
                if hasmutation_gp[i,j] == True :
                    if calls_gp[i,j] < 4:
                        mut_annotations['nts'] = mut_annotations['nts'] + nts[calls_gp[i,j]]
                        mut_annotations['nt_alt'] = nts[calls_gp[i,j]]
                        if len(mut_annotations['AA']) == 4:
                            mut_annotations['AA_gt'] = mut_annotations['AA_gt'] + mut_annotations['AA'][ maNT_gp[i,j] ]
                    elif calls_gp[i,j] == -1:
                        #if diverse (calls not a mutation), add minor and major call
                        mut_annotations['AA_gt'] = mut_annotations['AA_gt'] + mut_annotations['AA'][ maNT_gp[i,j] ]
                        mut_annotations['AA_gt'] = mut_annotations['AA_gt'] + mut_annotations['AA'][ minorNT_gp[i,j] ]
            if len(mut_annotations['AA']) == 4:
                mut_annotations['type'] = 'S' # eventually overwritten below if N
            # remove duplicates
            mut_annotations['nts'] = ''.join(OrderedDict.fromkeys( mut_annotations['nts'] ).keys()) # get only unique Nuc and keep order
            mut_annotations['AA_gt'] = ''.join(OrderedDict.fromkeys( mut_annotations['AA_gt'] ).keys()) # get only unique AA and keep order
            # record if nonsynonymous mutation
            if len(mut_annotations['AA_gt'])>1:
                mut_annotations['type'] = 'N'
            # Record all observed mutations across all isolates; E.g. K134Y, W47A, etc.
            if len(mut_annotations['AA_gt'])>1:
                mut_annotations['muts'] = []
                ancAA = mut_annotations['AA_gt'][0]
                derAAs = mut_annotations['AA_gt'][1:]
                for j,derAA in enumerate(derAAs):
                    mut_annotations['muts'].append( ancAA+str(mut_annotations['aa_pos'])+derAA )
            else:
                mut_annotations['muts'] = "."
            mut_annotations['NonSyn'] = len(mut_annotations['AA_gt'])>1
        else: #intergenic
            if int(mut_genenum[i])>0: # get info for gene prior SNP (if any)
                p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])-1] # both -1 necessary bcs list of df and df 0-based
                mut_annotations['gene1'] = p_anno.loc['gene']
                mut_annotations['locustag1'] = p_anno.loc['locustag']
                if 'orthologtag' in p_anno:
                    mut_annotations['orthologtag1'] = p_anno.loc['orthologtag']
                mut_annotations['product1'] = p_anno.loc['product']
                mut_annotations['distance1'] = mut_annotations['pos'] - p_anno.loc['loc2']
                if p_anno.loc['strand'] == -1:
                    mut_annotations['distance1'] = mut_annotations['distance1'] * -1

            if int(mut_genenum[i]+0.5) <= annotation_genes[p_chr-1].shape[0] and annotation_genes[p_chr-1].shape[1] !=0: # get info gene after SNP (if any); second conditional to evade empty chr
                p_anno = annotation_genes[p_chr-1].iloc[int(mut_genenum[i])] # -1 necessary bcs list of df 0-based; gene_id 0-based by we want following
                mut_annotations['gene2'] = p_anno.loc['gene']
                mut_annotations['locustag2'] = p_anno.loc['locustag']
                if 'orthologtag' in p_anno:
                    mut_annotations['orthologtag2'] = p_anno.loc['orthologtag']
                mut_annotations['product2'] = p_anno.loc['product']
                mut_annotations['distance2'] = p_anno.loc['loc1'] - mut_annotations['pos'] +1 # +1 to get correct bp distance
                if p_anno.loc['strand'] == 1:
                    mut_annotations['distance2'] = mut_annotations['distance2'] * -1
            # append_annotations intragenic
            #print(mut_annotations['distance1'])
            if ( 'distance1' in mut_annotations and mut_annotations['distance1'] > (-1*promotersize) and mut_annotations['distance1'] < 0) or ( 'distance2' in mut_annotations and mut_annotations['distance2'] > (-1*promotersize) and mut_annotations['distance2'] < 0):
                mut_annotations['type'] = 'P'
            else:
                mut_annotations['type'] = 'I'
            # define nts (repeat of intragenic code)
            if np.unique(refnti_gp[i,:])[0] != 4: # only if ref defined; [0] ok...see next line comment
                mut_annotations['nt_ref'] = nts[np.unique(refnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
                mut_annotations['nts'] = mut_annotations['nt_ref']
                if np.unique(ancnti_gp[i,:])[0] != 4:
                    mut_annotations['nt_anc'] = nts[np.unique(ancnti_gp[i,:])[0]] # refnti_gp should by definition never have > 1 unique element per row! >> [0] apropriate
            # extract derived genotype(s) across all samples
            for j,callidx in enumerate(calls_gp[i,:]):
                #print(str(j)+" "+str(callidx))
                #fixed mutation
                if hasmutation_gp[i,j] == True:
                    if calls_gp[i,j] < 4:
                        mut_annotations['nts'] = mut_annotations['nts'] + nts[calls_gp[i,j]]
                        mut_annotations['nt_alt'] = nts[calls_gp[i,j]]
            # remove duplicates
            mut_annotations['nts'] = ''.join(OrderedDict.fromkeys( mut_annotations['nts'] ).keys()) # get only unique Nuc and keep order
        lod_mutAnno.append(mut_annotations)
    dataframe_mut = pd.DataFrame(lod_mutAnno)
    return dataframe_mut

def annotate_sampleNames(samplenames,locations_long_names,patients_sbj,visits_sbj,locations_sbj):
    # extend sample name with patient/visit/location identifier. all in same order!
    extendend_sampleNames = np.copy(samplenames)
    for i,name in enumerate(extendend_sampleNames):
        extendend_sampleNames[i] = "S"+patients_sbj[i]+"_V"+str(visits_sbj[i])+"_"+locations_long_names[locations_sbj[i]]+"_"+name
    return extendend_sampleNames
       
def write_calls_sampleName_to_fasta(calls_for_tree,treeSampleNames,timestamp):
    fa_file = open(timestamp+".fa", "w")
    for i,name in enumerate(treeSampleNames):
        nucl_string = "".join(list(calls_for_tree[:,i]))
        fa_file.write(">" + name + "\n" + nucl_string + "\n")
    fa_file.close()

def generate_tree(calls_for_tree,treeSampleNamesLong,sampleNamesDnapars,refgenome,filetag,buildTree=False,writeDnaparsAlignment=False):
    # Write alignment file (as fasta)
    # calc NJ or Parsimonous tree or None
    # writeDnaparsAlignment==True for writing dnapars input for usage on cluster
    ts = time.time()
    timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
        
    # write alignment fasta,read alignment
    write_calls_sampleName_to_fasta(calls_for_tree,treeSampleNamesLong,timestamp) #timestamp.fa
    if writeDnaparsAlignment:
        # change tip labels and write phylip
        write_calls_sampleName_to_fasta(calls_for_tree,sampleNamesDnapars,timestamp+"_"+filetag+"_dnapars") #timestamp_dnapars.fa > for dnapars...deleted later
        # turn fa to phylip and delete fasta with short tip labels    
        aln = AlignIO.read(timestamp+"_"+filetag+"_dnapars.fa", 'fasta')
        AlignIO.write(aln, timestamp+"_"+filetag+".phylip", "phylip")
        subprocess.run(["rm -f " + timestamp+"_"+filetag+"_dnapars.fa"],shell=True)
        print("Written file: " + timestamp+"_"+filetag+".phylip")
        # write parameter file
        with open(timestamp+"_"+filetag+"_options.txt",'w') as file:
            file.write(timestamp+"_"+filetag+".phylip"+"\n")
            file.write("f"+"\n")
            file.write(timestamp+"_"+filetag+"_out.dnapars"+"\n")
            file.write("5"+"\n")
            file.write("V"+"\n")
            file.write("1"+"\n")
            file.write("y"+"\n")
            file.write("f"+"\n")
            file.write(timestamp+"_"+filetag+".tree"+"\n"+"\n")



    if buildTree=='PS':      
        # write phylip file with dnaparse compatible 10c samplenames
        write_calls_sampleName_to_fasta(calls_for_tree,sampleNamesDnapars,timestamp+"_dnapars") #timestamp_dnapars.fa > for dnapars...deleted later
        # turn fa to phylip and delete fasta with short tip labels    
        aln = AlignIO.read(timestamp+"_dnapars.fa", 'fasta')
        AlignIO.write(aln, timestamp+".phylip", "phylip")
        subprocess.run(["rm -f " + timestamp+"_dnapars.fa"],shell=True)
            
        # find dnapars executable
        dnapars_path = glob.glob('dnapars')
        path_extension = "../"
        backstop = 0
        while len(dnapars_path) == 0 and backstop <= 5:
            dnapars_path = glob.glob(path_extension+'dnapars')
            path_extension = path_extension + "../"
            backstop = backstop + 1
        if len(dnapars_path) == 0:
            raise ValueError('dnapars executable could not be located.')
        elif dnapars_path[0]=='dnapars':
            dnapars_path[0] = './dnapars'
        # write parameter file
        with open(timestamp+"_options.txt",'w') as file:
            file.write(timestamp+".phylip"+"\n")
            file.write("f"+"\n")
            file.write(timestamp+"_out.dnapars"+"\n")
            file.write("5"+"\n")
            file.write("V"+"\n")
            file.write("1"+"\n")
            file.write("y"+"\n")
            file.write("f"+"\n")
            file.write(timestamp+".tree"+"\n"+"\n")

        # run dnapars
        print("Build parsimony tree...")
        #print( dnapars_path[0] + " < " + timestamp+"_options.txt > " + timestamp+"_dnapars.log")
        subprocess.run([ "touch outtree"  ],shell=True)
        subprocess.run([ dnapars_path[0] + " < " + timestamp+"_options.txt > " + timestamp+"_dnapars.log"  ],shell=True)
        # print('done')
        # re-write tree with new long tip labels        
        tree = Phylo.read(timestamp+".tree", "newick")
        for leaf in tree.get_terminals():
            # print(leaf.name)
            idx = np.where(sampleNamesDnapars==leaf.name) 
            if len(idx[0]) > 1:
                warnings.warn("Warning: dnapars 10c limit leads to ambigous re-naming for "+leaf.name)
                idx = idx[0][0] #np.where returns: tuple with array with index
            else:
                idx = idx[0][0] #np.where returns: tuple with array with index
            leaf.name = treeSampleNamesLong[idx]
        Phylo.write(tree, timestamp+".tree", 'nexus')
        Phylo.write(tree, filetag+"_latest.nwk.tree", 'newick')
    
        # clean up
        subprocess.run(["rm -f " + timestamp+".phylip " + timestamp+"_options.txt " + timestamp+"_dnapars.log"],shell=True)

    elif buildTree == 'NJ':
        ## biopython tree build
        print("Build NJ tree...")
        # build starting tree (NJ)
        aln = AlignIO.read(timestamp+'.fa', 'fasta')
        calculator = DistanceCalculator('identity')
        constructor = DistanceTreeConstructor(calculator, 'nj')
        treeNJ = constructor.build_tree(aln)
        # Phylo.draw(treeNJ)
        Phylo.write(treeNJ,timestamp+"_NJ.tree","nexus")
        # build parsimonous tree
        #scorer = ParsimonyScorer()
        #searcher = NNITreeSearcher(scorer)
        #constructor = ParsimonyTreeConstructor(searcher, treeNJ)
        #treePS = constructor.build_tree(aln)
        #Phylo.write(treePS,timestamp+"_PS.tree","nexus")
    return timestamp

def build_table_for_tree_labeling(p_chr_table,treeSampleNamesLong,calls_for_tree,patient=""):
    # make new folder ('tree_counting'), wipe all previous content, add table
    if patient != "":
        subprocess.run(["rm -fr tree_counting/subject"+patient+" ; mkdir tree_counting/subject"+patient ],shell=True)
        with open("tree_counting/subject"+patient+"/for_tree_labeling.csv",'w') as csv_file:
            csv_file.write(",".join(np.append(np.array(['chr','pos']),treeSampleNamesLong))+"\n") # write header
            for i in range(p_chr_table.shape[0]):
                csv_file.write(",".join(np.append( np.array([str(p_chr_table[i,0]),str(p_chr_table[i,1])]) ,calls_for_tree[i,]))+"\n")
    else:
        subprocess.run(["rm -fr tree_counting/ ; mkdir tree_counting/ " ],shell=True)
    # build table    
    with open("tree_counting/for_tree_labeling.csv",'w') as csv_file:
        csv_file.write(",".join(np.append(np.array(['chr','pos']),treeSampleNamesLong))+"\n") # write header
        for i in range(p_chr_table.shape[0]):
            csv_file.write(",".join(np.append( np.array([str(p_chr_table[i,0]),str(p_chr_table[i,1])]) ,calls_for_tree[i,]))+"\n")


"""
TDL countMutations.py below until countMutations()
"""

def mutation_count(inputtree, lca, pos):
	emptycount=0
	ambiguous=['R','Y','M','K','S','W']
	ACGT=['A','C','T','G']

	tree=''
	for i,c in enumerate(inputtree):
		if c in ['(', ',', ')']:
			tree+=c
		elif c =='*':
			tree+=inputtree[i+1]
	
	#print tree
	#remove all question marks and ambiguities from the tree
	qfound=1
	while qfound==1:
		i=0
		qfound=0
		while i < len(tree):
			if tree[i]=='?' or tree[i] in ambiguous:
				qfound=1
				emptycount+=1
				left=i-1; balance=1; right=i+1
				if i==(len(tree)-2): 
					tree=tree[:-3]+')'
				elif i==1 and tree[i+3]==',':
					tree=tree[:i]+tree[i+2:]
				elif tree[i-1]==',' and tree[i+1]==',':
					tree=tree[:i-1]+tree[i+1:]
				#elif tree[i-1]=='(' and tree[i+1]==',' and tree[i+3]==',':
				#	tree=tree[:i]+tree[i+2:]
				elif i==0:
					tree=tree[2:]	
				elif tree[i+1]==')': #find open parenthesis side to remove
					while balance > 0:
						left+=-1
						if tree[left]=='(':
							balance+=-1
						if tree[left]==')':
							balance+=1
					tree=tree[:left]+tree[left+1:i-1]+tree[right+1:]
				elif tree[i-1]=='(': #find close parenthesis side to remove
					while balance > 0:
						right+=1
						if tree[right]=='(':
							balance+=1
						if tree[right]==')':
							balance+=-1
					if right == len(tree)-1:
						tree=tree[:left]+tree[i+2:-1]
					else:
						tree=tree[:left]+tree[i+2:right]+tree[right+1:]
			i=i+1
	
	#simplify tree while counting mutations. give preference to simplifying first.
	simplified=1
	mutlist={('A','C'):0,('G','C'):0,('T','C'):0,('A','G'):0,('C','G'):0,('T','G'):0, ('A','T'):0,('G','T'):0,('C','T'):0, ('T','A'):0,('G','A'):0,('C','A'):0}
	while simplified==1:
		#print tree
		#print mutlist
		i=0
		simplified=0
		while i < len(tree)-5:  #changed from i < len(tree)-4 1/25/2019 TDL -- may need to revert
			if tree[i]=='(' and tree[i+4]==')' and tree[i+1]==tree[i+3]: #we have pair of two closest related strains, simplify
				tree=tree[:i]+tree[i+3]+tree[i+5:]
				simplified=1
				i=i+1
			elif tree[i] in ACGT and tree[i]==tree[i+2]: #we have doublet of identical strains A,A, that aren't in a parenthesis-- reduce
				tree=tree[:i]+tree[i+2:]
				simplified=1
				i=i+1
			elif tree[i]=='(' and tree[i+2]==',' and tree[i+4]==',' and tree[i+1]==tree[i+5] and tree[i+1]!=tree[i+3] and tree[i+3] in ['A','C','T','G']: #we have (A,T,A ---> (A,A ... 
				anc=tree[i+1]
				mutlist[anc,tree[i+3]]+=1
				tree=tree[:i+3]+tree[i+5:]
				simplified=1
				i=i+10;
			elif tree[i]=='(' and tree[i+4]==')' and tree[i+1] in ACGT and tree[i+3] in ACGT: #we have a mutation, count it (A,G)
				if tree[i-2] in ACGT and tree[i+5]==')':
					anc=tree[i-2]
					if anc!=tree[i+1] and anc!=tree[i+3]: #rare, but has happened, when tree is... (A,(C,T))... and it is complicated to solve
						if tree[i-5] in ACGT:
							anc=tree[i-5]
						elif tree[i+7] in ACGT:
							anc=tree[i+7]
						if tree[i-2] != anc:
							mutlist[anc,tree[i-2]]+=1	
					if tree[i+1]!=anc:
						mutlist[anc,tree[i+1]]+=1
					if tree[i+3]!=anc:
						mutlist[anc,tree[i+3]]+=1		
					tree=tree[:i-2]+anc+tree[i+6:]
					simplified=1
				elif tree[i+6] in ACGT and tree [i-1]=='(':
					anc=tree[i+6]
					if anc!=tree[i+1] and anc!=tree[i+3]:  #rare, but has happened, when tree is... (A,(C,T)).
						if tree[i+9] in ACGT:
							anc=tree[i+9]
						elif tree[i-3] in ACGT:
							anc=tree[i-3]
						if tree[i+6] != anc:
							mutlist[anc,tree[i+6]]+=1
					if tree[i+1]!=anc:
						mutlist[anc,tree[i+1]]+=1
					if tree[i+3]!=anc:
						mutlist[anc,tree[i+3]]+=1		
					tree=tree[:i-1]+anc+tree[i+8:]
					simplified=1
					#print 'here3'
				elif tree[i-2] in ACGT and tree[i+6] in ACGT and tree[i-2]==tree[i+6] : #when A,(A,G),A
					anc=tree[i-2]
					if tree[i+1]!=anc:
						mutlist[anc,tree[i+1]]+=1
					if tree[i+3]!=anc:
						mutlist[anc,tree[i+3]]+=1		
					tree=tree[:i-2]+anc+tree[i+5:]
					simplified=1
				elif tree[i+6]=='(' and tree[i+10]==')' and tree[i+7]!=tree[i+9]: #this is a rare circumstance, I hope! When the tree is ...(A,T), (A,T)... count 2 mutations and remove both
					if tree[i+13] in ACGT:
						anc=tree[i+13]
						if tree[i+1]!=anc:
							mutlist[anc,tree[i+1]]+=1
						if tree[i+3]!=anc:
							mutlist[anc,tree[i+3]]+=1
						if tree[i+7]!=anc:
							mutlist[anc,tree[i+7]]+=1
						if tree[i+9]!=anc:
							mutlist[anc,tree[i+9]]+=1	
						tree=tree[:i-2]+anc+tree[i+15:]
						simplified=1	
						#print 'here3'
					elif tree[i-3] in ACGT:
						anc=tree[i-3]
						if tree[i+1]!=anc:
							mutlist[anc,tree[i+1]]+=1
						if tree[i+3]!=anc:
							mutlist[anc,tree[i+3]]+=1
						if tree[i+7]!=anc:
							mutlist[anc,tree[i+7]]+=1
						if tree[i+9]!=anc:
							mutlist[anc,tree[i+9]]+=1	
						tree=tree[:i-4]+anc+tree[i+13:]
						simplified=1
						#print 'here4'
					i=i+10 #try to simplify things to the both sides before dealing with something like ((A,G),((A,G),(A,C)))
					
			i+=1

#	print tree
	
	alreadyadded=[]
	#print(tree)
	for i,n in enumerate(tree):
		if n in ACGT and n!=lca and lca!='?' and n not in alreadyadded:
			mutlist[lca,n]+=1
			alreadyadded.append(i)
		elif i < (len(tree)-5) and tree[i]=='(' and tree[i+2]==',' and tree[i+4]==',' :
			if tree[i+1]==tree[i+5] and tree[i+1]!=tree[i+3] and tree[i+3] in ACGT: #we have (A,T,A, ... 
				anc=tree[i+1]
				mutlist[anc,tree[i+3]]+=1
				alreadyadded.append(i+3)


	
	
	muts=0
	for pair in mutlist.keys():
		muts+=mutlist[pair]
	
	#if muts==0:
		#print(tree)
	
	return muts, mutlist




def mutationtypes(tree, chart):


	ATCG=['A','C','T','G']
	f=open(chart,'r').readlines()

	
	
	for i, line in enumerate(f[1:]):
		#print(i)
		l=line.strip().split(',')
		if len(l) < 5:
			#print(l)
			continue
		chromosome=l[0]
		pos=l[1]
		
		#use first strain as lca
		lca=l[2]

				
		#count mutations
		newtree = annotateSNP(tree, chart, chromosome, pos)
		a, mutlist= mutation_count(newtree, lca, pos)
		#if a==0:
			#print('NO MUTS:')
			#print(chromosome, pos)
		#save trees
		if a != 1:
			f1=open(str(a)+'_'+chromosome+'_'+pos+'.tree','w')
		else:
			f1=open('1_'+chromosome+'_'+pos+'.tree','w')
		t=open('tempnexus.txt','r').readlines()
		for line in t:
			f1.write(line)
	

def annotateSNP(tree, dict, chromosome, pos):

	f=open(dict).readlines()
	fo=open('temptree.txt','w')
	
	header=f[0].strip().split(',')
	locations=[]
	
	for n,i in enumerate(header):
		if (i.startswith('S') or i.startswith('D') or i.startswith('R')) and len(i)<100:
			locations.append(n)
					
	for line in f:
		l=line.strip('\n').split(',')
		if l[0]==chromosome and l[1]==pos:
			for i in locations:
				if len(l) > i and len(l[i])>0:
					fo.write(header[i]+'\t'+l[i]+'\n')
				else:
					fo.write(header[i]+'\t?\n')
				#if i > len(l):
					#print(line, l, i)
			break

	fo.close()
	newtree = nameswap(tree,'temptree.txt')
	return newtree


def nameswap(tree, dictionary):

	f=open(dictionary).readlines()
	numStrains=len(f)
	dict={}
	annotation={}
	newname={}
	
	#print header for tempnexus
	fo=open('tempnexus.txt','w')
	fo.write('#NEXUS\nbegin taxa;\n\tdimensions ntax='+str(numStrains+1)+';\n\ttaxlabels\n')			
	colors={'A':'[&!color=#-16776961]', 'C':'[&!color=#-16725916]','G':'[&!color=#-3670016]','T':'[&!color=#-3618816]','?':'[&!color=#-16777216]'}
	ambiguous=['R','Y','M','K','S','W']
	for each in ambiguous:
		colors[each]='[&!color=#-16777216]'
				
	#get annotations
	f=open(dictionary, 'r').readlines()
	for line in f:
		if not line.startswith('#'):
			l=line.split()
			annotation[l[0]]=l[1]
	
	
	#combine names and annotations	
	for i in annotation.keys():
		newname[i]=i+ '--*'+ annotation[i] 
		#if i in dict.keys():
			#newname[i]=dict[i]+ '--*'+ annotation[i] #for dating newname[i]=dict[i]+ '--*'+ annotation[i]
		#else:
			#newname[i]= i + '--*'+ annotation[i] #for reference, etc.

	#print(newname)

	#make new tree
	f=open(tree,'r').readlines()
	
	newtree=''
	
	for line in f:
	#	print(line)
		i=0
		while i < len(line):
			#print(line[i])
			if line[i] not in ['S']: #S for to_mark_strain_name...
				newtree+=line[i]
				i+=1
			else:
				remainder=line[i:]
				nameend=remainder.find(':')
				oldname=remainder[:nameend]
				i=i+nameend
				if oldname in newname.keys():
					new=newname[oldname]
				else:
					new=oldname
				if new[-1]=='N':
					new=new[:-1]+'?'
				newtree+=new
				#print(newname)
				#print(oldname)
				fo.write('\t\''+new+'\''+colors[new[-1]]+'\n') #write down new colors to nexus file
	
	#write rest of nexus file
	fo.write('\t\'Reference\'[&!color=#-16777216]\n;\nend\n\nbegin trees;\n\ttree tree_1=[&R] ')
	for line in newtree:
		fo.write(line)
	fo.write('end;\n\nbegin figtree;\n')
	fo.write('\tset tipLabels.fontSize=10;\n')
	fo.write('end;')
	fo.close()
	return newtree	
# End TDL tree_counting functions

def countMutations(tree,chart):
    print('NOTE: Uses first strain as LCA')
    rc('font', **{'sans-serif':'Arial'})
    mutationtypes(tree,chart)
    #print('Done')

def plot_coverage_fwd_rev_stacked(chr_pos_gp,annotation_mutations,lod_fwd_cov,lod_rev_cov,timestamp,subject=""):
    # reconstruct TDL coverage plot for every SNP across all samples
    # pdfs: pdf/coverage_snp_fwd_rev/timestamp
    subprocess.run(["mkdir -p pdf/coverage_snp_fwd_rev/" + timestamp + "_"+subject+ " ; rm -f pdf/coverage_snp_fwd_rev/" + timestamp + "_"+subject + "/* "],shell=True)
    os.chdir('pdf/coverage_snp_fwd_rev/' + timestamp+ "_"+subject)
    for i in range(chr_pos_gp.shape[0]):
        # get pos,chr,locID,N/S/I/P, FQ for each positon and plot
        chr = chr_pos_gp[i,0]
        pos = chr_pos_gp[i,1] + 1 # turn 1-based bcs annotation_mutation is 1-based
        bool_pos_anno_df = (annotation_mutations['chr'] == chr) &  (annotation_mutations['pos'] == pos)
        if annotation_mutations.loc[ bool_pos_anno_df , 'locustag' ].isnull().any(axis=0):
            locID = "nan"
        else:
            locID = annotation_mutations.loc[ bool_pos_anno_df , 'locustag' ].values[0] # locustag
        anno = annotation_mutations.loc[ bool_pos_anno_df , 'type' ].values[0] # N/S/P/I
        qual = str(int(annotation_mutations.loc[ bool_pos_anno_df , 'quals' ].values[0])) # single qual value (based on FQ: best of lowest scoring pair)
        identifier_string = str(chr) + "_" + str(pos) + "_" + locID + "_" + anno + "_" + qual        
        print("Plot coverage fwd/rev: "+identifier_string )
        plot_stacked_paired( [lod_fwd_cov[i],lod_rev_cov[i]] , labels=["fwd", "rev"] , title=identifier_string )
    os.chdir('../../../')

def build_dataframe_coverage_info(goodpos2useTree,NTs,SampleNames,maNT,minorNT,coverage_forward_strand,coverage_reverse_strand,maf,minorAF):
    listOfDF_fwd = []
    listOfDF_rev = []
    NTs_wN = np.append(NTs,'N') # there are N in maf/ninorNT
    for i in goodpos2useTree: # array contains index of p
        lod_fwd = []
        lod_rev = []
        # extract for each pos-sample pair the read_count for ACTG (we record maf & minor >> thus max 2 nucl per pos)
        for j,sample in enumerate(SampleNames):
            spl_fwd = {}
            spl_rev = {}
            spl_fwd['sample'] = sample
            spl_rev['sample'] = sample
            for k,nuc in enumerate(NTs_wN):
                # dict had to be defined within if else, otherwise get overwritten
                if k == maNT[i,j]:
                    # print('maNT',k,nuc)
                    spl_fwd[nuc] = ( coverage_forward_strand[i,j] * maf[i,j] ) # recalc #Nucleotide calls based on coverage and major allele freq
                    spl_rev[nuc] = ( coverage_reverse_strand[i,j] * maf[i,j] ) # recalc #Nucleotide calls based on coverage and major allele freq
                elif k == minorNT[i,j]:
                    # only report minorAF when majorNT not set 4 during filtering!
                    if maNT[i,j] != 4:
                        spl_fwd[nuc] = ( coverage_forward_strand[i,j] * minorAF[i,j] ) # recalc #Nucleotide calls based on coverage and major allele freq
                        spl_rev[nuc] = ( coverage_reverse_strand[i,j] * minorAF[i,j] ) # recalc #Nucleotide calls based on coverage and major allele freq
                else:
                    spl_fwd[nuc] = 0
                    spl_rev[nuc] = 0
            lod_fwd.append(spl_fwd)
            lod_rev.append(spl_rev)
            df_fwd = pd.DataFrame(lod_fwd)
            df_fwd = df_fwd.set_index('sample')
            df_rev = pd.DataFrame(lod_rev)
            df_rev = df_rev.set_index('sample')
        listOfDF_fwd.append(df_fwd)
        listOfDF_rev.append(df_rev)
    return [listOfDF_fwd,listOfDF_rev]

def plot_stacked_paired(dfall, labels=None, title="notitle",  H="/", **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
labels is a list of the names of the dataframe, used for the legend
title is a string for the title of the plot
H is the hatch used for identification of the different dataframe"""
    #NOTE: when major NT == 'N', minor allele frequency not shown. Column only represents major AF in this specific case. Otherwise column presents majorAF plus minorAF
    n_df = len(dfall)
    n_col = len(dfall[0].columns) 
    n_ind = len(dfall[0].index)
    
    # define plot width, min = 7
    plot_width = int(n_ind/5)
    if plot_width < 7:
        plot_width = 7
    fig = plt.figure(figsize=(plot_width, 8))
    axe = plt.subplot(111)

    for df in dfall : # for each data frame
        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ylim=(0,100),
                      ax=axe,
                      legend=False,
                      grid=False,
                      **kwargs)  # make bar plots

    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                rect.set_hatch(H * int(i / n_col)) #edited part     
                rect.set_width(1 / float(n_df + 1))

    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    axe.set_xticklabels(df.index, rotation = 90)
    axe.set_title(title)

    # Add invisible data to add another legend
    n=[]        
    for i in range(n_df):
        n.append(axe.bar(0, 0, color="gray", hatch=H * i))

    l1 = axe.legend(h[:n_col], l[:n_col], loc=[1.01, 0.5])
    if labels is not None:
        l2 = plt.legend(n, labels, loc=[1.01, 0.1]) 
    axe.add_artist(l1)
    fig.savefig(title + ".pdf")
    plt.close()

    
def parallel_evolution_counting_and_simulation(num_mutations_genome, num_mutations_genic , mutation_number_threshold , mutation_density_threshold , numtrials, max_muts_per_gene_to_track, chr_pos_gp, ref_genome_folder, annotation_genes):
    ## Parallel evolution counting and simulation
    # Inputs: genome, number of expected mutations in the genome
    # Output: number of genes with a mutation density (mutations per length of
    # gene) over some threshold

    # num_mutations_genome=number_of_mutations_on_genome
    # num_mutations_genic=number_of_mutations_genic
    # mutation_number_threshold=params_dict['Min_num_mutations_cand']
    # mutation_density_threshold=params_dict['Min_mutation_density_cand']
    # numtrials=params_dict['NumTrialsSim']
    # chr_pos_gp=chr_pos_gp
    # ref_genome_folder=params_dict['ref_genome_folder']
    


    # Get data structures for analysis
    [chrstarts, genomelength,scafnames]= genomestats(ref_genome_folder);
    [locus_tag,chr_tag] = genomic_position_all(annotation_genes, genomelength, chrstarts) #loc_tag numbers all genes(CDS) across genome (intergenic/rna:0.5); chr_tag: numbers genes per chr (intergenic: #chr+0.5)
    #    locus_tag >> cds_number_by_position_tags
    #    chr_tag >> cds_number_by_position_indices
    # extract vector genelengths that follows order
    
    annotation_genes_sglDF = pd.concat(annotation_genes,sort=False)
    start_pos_gene = annotation_genes_sglDF['loc1'].values
    start_pos_gene = start_pos_gene[~np.isnan(start_pos_gene)] # remove nan, caused by empty chr
    end_pos_gene = annotation_genes_sglDF['loc2'].values
    end_pos_gene = end_pos_gene[~np.isnan(end_pos_gene)] # remove nan, caused by empty chr
    # gene lengths. incl rna genes (however when used those indexes not present in locus_tag (0.5))
    genelengths = (end_pos_gene - start_pos_gene)
    
    expectedNumberGenesMultipleMutations = np.zeros(numtrials) #initialize vector to store simulation results
    expectedNumberOfGenesWithNmutations = np.zeros((numtrials,max_muts_per_gene_to_track));

    # start sim
    for i in range(numtrials):
    
        # Pick random positions on the genome to mutate
        # Initially get 10x the number of positions you need --> later filter only for num_mutations_genic mutations
        randpos = randint(1, genomelength, 10*num_mutations_genome , dtype=int)
        # Does NOT assume any mutational spectrum!!!
        
        #TDL has scripts for doing the same thing for operon and pathway
        #levels, but the inputs they are generated on may not be relevant for
        #your genome
    
        # Find out in which genes these mutations occurred
        genenums = locus_tag[randpos];
        genenums = genenums[genenums>0.5]; #remove intergenic/rna-gene mutations
        genenums = genenums[0:num_mutations_genic]; # only take as many as you need (generated extra because some removed due to not being on a gene)
        # Get a list of the genes mutated, along with the number of times each
        [sim_unique_genes , sim_mutspergene] = np.unique(genenums,return_counts=True)

    
        # This calculates the number of mutations per gene length for each gene
        # on which there were simulated mutations
        sim_mutspergenelength = sim_mutspergene/genelengths[sim_unique_genes.astype(int)-1]; #-1 bcs sim_unique_genes 1-based tagnumber > 1st gene is position 0 in genelengths

        # The final step finds the number of genes that were mutated multiple 
        # times and above the threshold mutation density (per unit length). 
        # Result saved, indexed by trial number.
        expectedNumberGenesMultipleMutations[i] = sum( (sim_mutspergenelength >= mutation_density_threshold) & (sim_mutspergene >= mutation_number_threshold)); # Arolyn, 2018.12.14 > to >= to match input from function
    
        # The second piece of information this script returns is the number of
        # genes with >m mutations
        for j in range(max_muts_per_gene_to_track):
            expectedNumberOfGenesWithNmutations[i,j] = sum( sim_mutspergene >= j+1 ) # +1 due to 0-based
            
    return [expectedNumberGenesMultipleMutations, expectedNumberOfGenesWithNmutations ]

def codon_composition_table( allbases, allcodons ):
    ''' Build table containing base counts (col) for each codon (row) '''
    codoncompositiontable = np.zeros((len(allcodons),len(allbases)),dtype=int)
    for i,codon in enumerate(allcodons):
        for j,base in enumerate(allbases):
            codoncompositiontable[i,j] = codon.count(base)
    return codoncompositiontable

def codon_mutation_table(allmuts , allcodons , codon_all_dict):
    ''' Generates table of probabilities that a given mutation on a codon is nonsynonymous '''    
    table = np.zeros( (len(allcodons), len(allmuts) ) ); # Initialize table
    for i,codon in enumerate(allcodons):
        for j,mut in enumerate(allmuts):
            # Calculates the probability that a given mutation is nonsynonymous on a given codon and then updates the table
            table[i,j] = prob_nonsyn_codon_mutation( codon, mut , codon_all_dict);
    return table

def prob_nonsyn_codon_mutation(codon,mut,codon_all_dict):
    ''' Calculates the probability that a given mutation leads to a nonsynonymous change across a given codon '''
    aa0 = codon_all_dict[codon] # AA of codon
            
    # Find the positions on the codon at which mutation could occur
    possiblemuts=[i for (i, base) in enumerate(codon) if base == mut[0]]
    ctr = 0
    if len(possiblemuts) == 0: # if the mutation cannot occur on this codon
        probability = float('nan');
    else: # mut can occur at least once
        for pos in possiblemuts:
            newcodon = list(codon)
            newcodon[pos] = mut[1] # mutate codon position that carries mut[0]
            aa1 = codon_all_dict[ "".join(newcodon) ]
            if aa0 != aa1:
                ctr += 1
        probability = ctr/len(possiblemuts) # fraction of mutations that were nonsynonymous
    return probability
     
def codons_in_genome(annotation_genes,allcodons):
    ''' Get probability for each codon across all CDS annotated in the reference genome (annotation_genes)'''
    # possibility to add a flag in order to restrict analysis to genomic region
    annotation_genes_sglDF = pd.concat(annotation_genes,sort=False)
    annotation_genes_sglDF_CDS = annotation_genes_sglDF.loc[ (annotation_genes_sglDF['type'] ==  'CDS') | (annotation_genes_sglDF['type'] ==  'gene') ]
    # Tally codon occurrences over all proteins in genome
    codonCounts = np.zeros( len(allcodons) ,dtype=int ); # for storing tally of codons
    for i, row in annotation_genes_sglDF_CDS.iterrows():
        seq = str(row['sequence'])
        startCodon = row['codon_start']
        startPos = startCodon * 3
        codons_gene = [seq[i:i+3] for i in range(startPos, len(seq), 3)]
        codons_gene = collections.Counter(codons_gene) # builds dictionary with all codons (key) in list and counts (value)
        for i,codon in enumerate(allcodons):
            if codon in codons_gene.keys():
                codonCounts[i] += codons_gene[codon]
    return codonCounts/sum(codonCounts) # probabilities of codons sorted by allcodons order
 
def mutation_probability(params_dict,annotation_genes):
    ''' This script examines the probability of a nonsynonymous mutation on a 
        reference genome given some mutational spectrum. '''
    # based on arolyn's m script: @Feb 2019
    ## Define DNA bases, possible mutations, and codons
    # All possible mutations to be considered:
    allbases = np.array(['A','T','G','C']) #fourDNAbases
    allmuts = np.array(['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG']); # 'AT' denotes A to T
    
    # All standard codons
    standard_codon_table = CodonTable.unambiguous_dna_by_id[1]
    allcodons = np.array([c for c in standard_codon_table.forward_table.keys()],dtype=object) # standard codon table, but miss stop codosn
    allcodons = np.append(allcodons , np.array([c for c in standard_codon_table.stop_codons],dtype=object) ) # 64 codons
    # build a combined dictionary (codon>>AA) containing regular and stop codon AA (stop:*)
    codon_all_dict = {}
    for c in allcodons:
        if c in standard_codon_table.forward_table.keys():
            codon_all_dict[c] = standard_codon_table.forward_table[c]
        else:
            codon_all_dict[c] = "*"

    # Generate table of codon composition by base
    codoncompositiontable = codon_composition_table( allbases, allcodons );
    # Rows (64) = all possible codons
    # Columns (4) = number of A's/T's/G's/C's in each codon
    
    ## Generate table of probabilities of nonsynonymous mutations
    # Rows (64) = all possible codons
    # Columns (12) = all possible mutations
    # Entries = probability of nonsynonymous mutation (given a mutation and a codon)
    # Note: If the given mutation cannot be applied to a given codon, the entry is nan. 

    codonnonsyntable = codon_mutation_table( allmuts, allcodons , codon_all_dict );

    ## Calculate mutational spectrum
    # Assuming all types of mutations are equally likely to occur:
    # Get mutational spectrum from experiments, but ignore if fwd or rev strand
    # allmuts = ['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG']
    #AT, TA  0
    #AC, TG  1
    #AG, TC  2
    #GC, CG  3
    #GT, CA  4
    #GA, CT  5
    if params_dict['substitution_spectrum'] == None:
        mutationalspectrum = [1/12] * 12 # uniform distribution. replace with time stratified observation based on all samples
    else:
        afile = open(params_dict['substitution_spectrum'], 'rb')
        mutationalspectrum = pickle.load(afile)
        afile.close()
        print('Observed substitution spectrum loaded')

    ## Calculate codon distribution in reference genome ordered as in allcodons
    codondistribution = codons_in_genome( annotation_genes , allcodons );
    
    ## Calculate probability of nonsynonymous mutation
    # Takes into account: mutation spectrum, abundance of codons on genome,
    # abundance of bases in codons, probability of a given mutation being
    # nonsynonymous on a givne codon...
    probnonsyn = 0; # probability of nonsynonymous mutations over all possible mutations

    for i,mut in enumerate(allmuts): # loop through all possible mutations
        # Probability that this mutation occurs
        prob_base_base = mutationalspectrum[i]
        # Find the codons that can undergo this mutation:
        base = mut[0] # base that gets mutated; ex. A
        baseindex = np.where(allbases==base) # ex. A is indexed in position 1
        basecodonoccurrences = codoncompositiontable[:,baseindex].flatten()  # ex. how many A's in each codon
        # Indices of codons that have the relevant initial base:
        basecodonoccurrences_bool = basecodonoccurrences > 0; # ex. AAT has A's but GGC does not
        # Probability that this mutation occurs on a given relevant codon
        # Take into account base composition of codons
        basecountincodon = basecodonoccurrences[ basecodonoccurrences_bool ];
        # Take into account codon abundance on reference genome
        probcodonongenome = codondistribution[ basecodonoccurrences_bool ]
        # Combine these two probabilities
        probmutoncodon = basecountincodon*probcodonongenome
        # Renormalize (sum = 1 over all relevant codons)
        probmutoncodon = probmutoncodon/sum(probmutoncodon);

        # Probability that this mutation is nonsynonymous at each relevant codon
        thismutnonsynoncodon = codonnonsyntable[:,i];
        probmutnonsynoncodon = thismutnonsynoncodon[basecodonoccurrences_bool];
        # Overall probability that this mutation is nonsynonymous over all possible codons
        probmutnonsyn=prob_base_base*sum(probmutoncodon*probmutnonsynoncodon);

        # Add contribution of this mutation to the total probability of a nonsynonymous mutation
        probnonsyn += probmutnonsyn  
        
    print('Probability of nonsynonymous mutation across genome: ' + str(probnonsyn) )
    return probnonsyn # Prob for N occuring 

def parallel_evo_module(goodpos2use,contig_positions,annotation_mutations, annotation_genes, params_dict,plot=True):
    ''' Module to calculate parallel evolution 
    # Test excess of genes with multiple mutations (thresholds for candidates defined in parameters dict )
    # Test excess of NonSyn (dNdS) in candidate genes compared to expectation given reference genome 
    # NOTE: dNdS uses simple substitution model (equal P per mut), which should be updated based on observations '''
    print('Parallel evolution inference.')
    # Find mutations that are adjacent to each other.
    # True for any mutation that is followed by an adjacent mutation, False if not
    # keep trailing bases of adjacent sets (2+)
    chr_pos_gp = contig_positions[goodpos2use,]
    bool_adjacent_mut = np.full( chr_pos_gp.shape[0] , False, dtype=bool) #mutated_genes_pos = []; # keep track of the positions of each event in the table
    for i in range(chr_pos_gp.shape[0]):
        chr = chr_pos_gp[i,0]
        pos = chr_pos_gp[i,1]
        if (i+1) <= (chr_pos_gp.shape[0]-1): # (chr_pos_gp.shape[0]-1) bcs shape[0] not 0-based
            if chr == chr_pos_gp[i+1,0] and (pos+1 == chr_pos_gp[i+1,1]):
                bool_adjacent_mut[i] = True
        
    # get info for candidate genes
    mutated_genes = np.zeros(0,dtype=int); # keep track of gene number
    mutated_genes_tally = np.zeros(0,dtype=int); # keep track of how many SNPs there are on this gene
    mutated_genes_lengths = np.zeros(0,dtype=int); # keep track of the length of this gene
    locustags_all = np.zeros(0,dtype=object) # record the locustag
    orthologtags_all = np.zeros(0,dtype=object) # record the orthologtags
    for i, row in annotation_mutations.iterrows():
        if bool_adjacent_mut[i] == False: # ignore leading adjacent mutations
            gene_num_global = row['gene_num_global']
            if gene_num_global != 0.5: # non-genic 0.5
                if gene_num_global in mutated_genes:
                    mutated_genes_tally[-1] = mutated_genes_tally[-1] + 1
                else:
                    mutated_genes = np.append( mutated_genes , gene_num_global )
                    mutated_genes_tally = np.append( mutated_genes_tally , 1 )
                    mutated_genes_lengths = np.append( mutated_genes_lengths , (row['loc2']-row['loc1']+1) ) # +1 bcs loc1/2 1-based, thus underestimate L by 1bp
                    locustags_all = np.append(locustags_all , row['locustag'])
                    orthologtags_all = np.append(orthologtags_all , row['orthologtag'])
    mutated_genes_tally_perGeneLen = mutated_genes_tally/mutated_genes_lengths;
        
    #% Define candidates for selection
    mutation_number_threshold = params_dict['Min_num_mutations_cand']; # minimum number of mutations per gene
    mutation_density_threshold = params_dict['Min_mutation_density_cand']; # minimum number of mutations per 1000 bp
    
    mutated_genes_of_interest = ( mutated_genes_tally >= mutation_number_threshold) & (mutated_genes_tally_perGeneLen >= mutation_density_threshold );
    num_mutated_genes_of_interest = sum(mutated_genes_of_interest);
    print('Number of genes with multiple mutations: ' + str(num_mutated_genes_of_interest))
    with open("parallel_evo_module_results" + params_dict['subjectID'] + ".txt", "w") as myfile:
        myfile.write('Number of genes with multiple mutations: ' + str(num_mutated_genes_of_interest) + "\n")
        myfile.write('Minimum number mutations required: ' + str(mutation_number_threshold) + "\n")
        myfile.write('Minimum SNP density for candidate genes: ' + str(mutation_density_threshold) + "\n")        
    
    # break if no candidate genes present
    if num_mutated_genes_of_interest == 0:
        print('No genes with multiple mutation found! >> skip adaptive evolution analysis')
        return [np.array([]),pd.DataFrame()] # return empty array and dataframe
    
    # get annotation_mutation for SNPs with signature of parallel evolution
    annotation_mutation_paraSignal = annotation_mutations.loc[ annotation_mutations['gene_num_global'].isin( mutated_genes[mutated_genes_of_interest] ) ]



    
    # =============================================================================
    # Test candidate genes for parallel mutation enrichment
    # =============================================================================
    # NOTE: Simulation based on random genic SNPs in genome and does not use a specific substitution model
    number_of_mutations_on_genome = len(goodpos2use)
    number_of_mutations_genic = sum(mutated_genes_tally); # only genic/CDS!
    
    [expectedNumberGenesMultipleMutations, expectedNumberOfGenesWithNmutations] = parallel_evolution_counting_and_simulation(number_of_mutations_on_genome,number_of_mutations_genic,params_dict['Min_num_mutations_cand'],params_dict['Min_mutation_density_cand'],params_dict['NumTrialsSim'],params_dict['max_muts_per_gene_to_track'],chr_pos_gp , params_dict['ref_genome_folder'],annotation_genes)
    # expectedNumberGenesMultipleMutations: num candidate genes per sim for parallel evolution based on params_dict['Min_num_mutations_cand'],params_dict['Min_mutation_density_cand']
    # expectedNumberOfGenesWithNmutations: row==sim; col==genes with mutations. col[0] == 1 mutation!, col[1]==2mut...
    
    simProbForObsCand = 1-sum(expectedNumberGenesMultipleMutations < num_mutated_genes_of_interest)/len(expectedNumberGenesMultipleMutations)
    print('Simulation-based probability to observe ' + str(num_mutated_genes_of_interest) + ' canddiate genes for parallel evolution is: '+str(simProbForObsCand))
    print('Simulated a mean number of candidate genes is ' + str(np.mean(expectedNumberGenesMultipleMutations)))
    with open("parallel_evo_module_results" + params_dict['subjectID'] + ".txt", "a") as myfile:
        myfile.write('Simulation-based probability to observe ' + str(num_mutated_genes_of_interest) + ' canddiate genes for parallel evolution is: '+str(simProbForObsCand) + "\n")
        myfile.write('Simulated a mean number of candidate genes is ' + str(np.mean(expectedNumberGenesMultipleMutations)) + "\n")

    # calc prob to observe candidate mut counts
    locustags_cand = locustags_all[mutated_genes_of_interest]
    orthologtags_cand = orthologtags_all[mutated_genes_of_interest]
    mut_cand_tally = mutated_genes_tally[mutated_genes_of_interest]
    prob_cand_nummut = np.ones(mut_cand_tally.shape)
    for i,nummut in enumerate(mut_cand_tally):
        prob_cand_nummut[i] = 1- (sum(expectedNumberOfGenesWithNmutations[:,0:(nummut-1)])/sum(expectedNumberOfGenesWithNmutations))
    res_cand_nummut = np.vstack((orthologtags_cand,locustags_cand,mut_cand_tally,prob_cand_nummut)).T
    print('Orthologtag Locustag NumMuts ProbSim')
    print(res_cand_nummut)    
    with open("parallel_evo_module_results" + params_dict['subjectID'] + ".txt", "a") as myfile:
        myfile.write('Orthologtag Locustag NumMuts ProbSim'  + "\n")
        np.savetxt(myfile,  res_cand_nummut,fmt="%s")

    # =============================================================================
    # dN/dS calculation
    # =============================================================================
    # calc if we observe more NS than expected (simulated)
    
    if np.mean(expectedNumberGenesMultipleMutations) <= num_mutated_genes_of_interest:
        print('Investigate dN/dS')
        # Observed dNdS
        ## Gene numbers of genes of interest
        cand_genes_geneNumGlobal = mutated_genes[mutated_genes_of_interest]; # these_gene_nums
        cand_muts_N = 0;
        cand_muts_S = 0;
        for global_tag in cand_genes_geneNumGlobal:
            cand_mut_anno = annotation_mutations.loc[ annotation_mutations['gene_num_global'] == global_tag ]
            for i,row in cand_mut_anno.iterrows():
                if row['type'] == "N":
                    cand_muts_N += 1
                elif row['type'] == "S":
                    cand_muts_S += 1
        
        if (cand_muts_N+cand_muts_S) == 0: # necessary in case no S/N mutations but all P/U
            print("No dN/dS calculated >> no S/N mutations.")
        else:
            probN_obs = cand_muts_N/(cand_muts_N+cand_muts_S)        
            # Simulate genome specific expected probability for NonSyn to occur
            # simulate only if value not yet calculated. > stored in regenome/probNsim.pk     
            # probN_sim = 0.6829; 
            probN_sim = mutation_probability(params_dict,annotation_genes)
            # NOTE: mutation_probability uses a pre-calculted substitution spectrum (mutationalspectrum) when specified in params or a uniform mutation spectrum.
            
            # Binomial test
            prob_obsN_excess = stats.binom_test(cand_muts_N, n=(cand_muts_N+cand_muts_S), p=probN_sim, alternative='greater')
            print('Observed P(N): ' + str(probN_obs) + '; N: ' + str(cand_muts_N) + ', S: ' + str(cand_muts_S))
            print('Expected P(N): ' + str(probN_sim) )
            print('Probability of enrichment in observed count of N: ' + str(prob_obsN_excess) ) 
            with open("parallel_evo_module_results" + params_dict['subjectID'] + ".txt", "a") as myfile:
                myfile.write('Observed P(N): ' + str(probN_obs) + '; N: ' + str(cand_muts_N) + ', S: ' + str(cand_muts_S) + "\n")
                myfile.write('Expected P(N): ' + str(probN_sim) + "\n")
                myfile.write('Probability of enrichment in observed count of N: ' + str(prob_obsN_excess) + "\n")
        

    # =============================================================================
    #     Plot
    # =============================================================================
    if plot:
        plt.rcParams.update({'font.size': 14}) # all label size
        f = figure()
        ax=f.add_subplot(111)
    
        # histogram depicting simulated distribution of expected candidate genes, filter criteria for candidates, P for observation cand count
        mybin = np.arange(max(np.unique(expectedNumberGenesMultipleMutations))+2)-0.5 # +2 in order to now loose last bin due to arange; -0.5 needed for bars align at center with x axis ticks
        plt.hist(expectedNumberGenesMultipleMutations,bins=mybin,rwidth=0.8,color='#607c8e')
        xticks(np.array([0,1,2]))
        plt.ylabel('Simulated counts')
        plt.xlabel('Number of genes with mutliple mutations')
        plt.axvline(x=len(mut_cand_tally),color='violet') 
        text(0.98, 0.95, "min #mut:"+str(params_dict['Min_num_mutations_cand'])+"; min density:"+str(params_dict['Min_mutation_density_cand']), fontsize=12,horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
        text(0.98, 0.88, "P("+ str(len(mut_cand_tally)) + ") = " + str(np.around(simProbForObsCand,3)), fontsize=12,horizontalalignment='right',verticalalignment='center',transform = ax.transAxes)
        #save pdf
        subprocess.run(["mkdir -p pdf/adaptive_evo/ "],shell=True)
        f.savefig('pdf/adaptive_evo/' + params_dict['timestamp'] + "_" + params_dict['subjectID'] + ".pdf")
        plt.close()
        print('Plotted: pdf/adaptive_evo/' + params_dict['timestamp'] + "_" + params_dict['subjectID'] + ".pdf")
    return [res_cand_nummut,annotation_mutation_paraSignal]


def get_index_visit_sampleNames(spl_names_long,treesampleNamesLong=True):
    # get idx V1,2,3,4,5 based on treesampleNamesLong
    # later extended funtion to also use visits variable instead of treesampleNamesLong
    list_visit_idx = []
    if treesampleNamesLong:
        for v in ['V1','V2','V3','V4','V5']:
            match_visit = np.array([i for i,x in enumerate(spl_names_long) if re.search('_'+v+'_',x)])
            list_visit_idx.append(match_visit)
    else:
        for v in [1,2,3,4,5]:
            match_visit = np.array([i for i,x in enumerate(spl_names_long) if x==v])
            list_visit_idx.append(match_visit)        
    return list_visit_idx

def calc_regression(mean_count_per_visit,subject_fld_label,num_samples,num_visits):
    # calculate regression 
    # return dictionary with results and some additional metainfo
    from scipy import stats
    import scipy
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(mean_count_per_visit[:,0] , mean_count_per_visit[:,1] )
    regress_dict = {}
    regress_dict['slope'] = slope
    regress_dict['intercept'] = intercept
    regress_dict['r_sq'] = r_value**2
    regress_dict['p_val'] = p_value
    regress_dict['std_err'] = std_err
    regress_dict['subject'] = subject_fld_label
    regress_dict['num_samples'] = num_samples
    regress_dict['num_visits'] = num_visits
    return regress_dict

def plot_molecular_clock(input_data,hdr,mean_count_per_visit,regress_dict,basescale,subject_identifier):
    ''' plot molecular clock for genome-size corrected inference. rescale to genome-wide clock values using eg. 'basescale'=100000 '''
    # turn to df for seaborn
    # basescale 
    data_df = pd.DataFrame(data=input_data,    # values
                 index=np.arange(input_data.shape[0]),    # 1st column as index
                 columns=hdr)  # 1st row as the column names        
    data_df['is_even'] = (data_df['Month'] % 2) == 0 # add T/F for coloring light/darkgrey
    xpos = np.unique(data_df['Month']) # time in months for adding line
    slope = regress_dict['slope']
    r_value_sq = regress_dict['r_sq']
    subject_id = regress_dict['subject']
    intercept = regress_dict['intercept']
    
    ## plot
    # !!! someone should add jitter!!!
    useylim = max(data_df['rate']*basescale)+max(data_df['rate']*basescale)*0.2
    
    plt.figure(figsize=(5,5))
    plt.scatter(data_df['Month'],data_df['rate']*basescale,facecolors='none', edgecolors='grey',alpha=.5)
    plt.scatter(mean_count_per_visit[:,0], mean_count_per_visit[:,1]*basescale, marker="x", s=100, c='black')
    plt.plot(xpos,(intercept*basescale) + slope*(xpos)*basescale,'b-')
    plt.ylim((0, useylim))
    plt.suptitle(subject_identifier+'\nMolecular. rate (mut/y/'+str(basescale)+'b) = '+str( np.round(12*slope*basescale ,3) ) + '\n r^2 = '+str( np.round(r_value_sq ,3) ), fontsize=12)
    plt.xlabel("Time (months)",fontsize=14)
    plt.ylabel("Molecular rate per "+str(basescale)+'b/y',fontsize=14)    
    subprocess.run(["mkdir -p pdf/molclock/ "],shell=True) # build output folder
    plt.savefig('pdf/molclock/' + subject_id + "_assembleBased_v2.pdf")
    plt.close()
    print('Molecular clock analysis: slope: ' + str(np.round(12*slope*basescale ,3)) + ", r^2: " + str( np.round(r_value_sq ,3) ) )
    print("Done. pdf/molclock/" + subject_id + "_assembleBased.pdf")



# def plot_molecular_clock(input_data,hdr,mean_count_per_visit,regress_dict,basescale = 1000000):
#     ## old version. deprecated. kept for legacy purpose only. if not your legancy please remove from analysispy_module.py
#     # turn to df for seaborn
#     data_df = pd.DataFrame(data=input_data,    # values
#                  index=np.arange(input_data.shape[0]),    # 1st column as index
#                  columns=hdr)  # 1st row as the column names        
#     data_df['is_even'] = (data_df['Month'] % 2) == 0 # add T/F for coloring light/darkgrey
#     xpos = np.unique(data_df['Month']) # time in months for adding line
#     slope = regress_dict['slope']
#     r_value_sq = regress_dict['r_sq']
#     subject_id = regress_dict['subject']
#     intercept = regress_dict['intercept']
    
#     ## plot    
#     sns.set(font_scale=1.3,style="white",rc={"lines.linewidth": 2.7}) 
#     ax = sns.catplot(x="Month", 
#                 y="rate", 
#                 data=data_df,
#                 kind="swarm",
#                 aspect=1.8, 
#                 hue='is_even',
#                 legend=False,
#                 palette=sns.color_palette(['gray', 'silver']),
#                 zorder=0)    
#     ax.set(ylim=(0, max(data_df['rate'])+max(data_df['rate'])*0.2),xlim=(-1,np.max(data_df['Month'])+2),ylabel="Number of Mutations",xlabel="Time (in months)")
#     ax = scatter(mean_count_per_visit[:,0], mean_count_per_visit[:,1], marker="x", s=100, c='black')
#     ax = sns.lineplot(xpos,
#                       intercept + slope*(xpos), 
#                       legend=False,
#                       color="black")
#     text(0.5, 0.04, ("R^2 = "+str( np.round(r_value_sq ,3) )), fontsize=15,horizontalalignment='left',verticalalignment='center',transform = ax.transAxes)
#     text(0.5, 0.12, ("Molecular rate (mut/year) = "+str( np.round(12*slope ,3) )), fontsize=17,horizontalalignment='left',verticalalignment='center',transform = ax.transAxes)
#     ax.set_title('Subject '+subject_id)
#     #save pdf
#     subprocess.run(["mkdir -p pdf/molclock/ "],shell=True)
#     fig = ax.get_figure()
#     plt.show()
    
#     fig.savefig('pdf/molclock/' + subject_id + "_assembleBased.pdf")
#     plt.close()
#     print('Molecular clock analysis: slope: ' + str(np.round(12*slope ,3)) + ", r^2: " + str( np.round(r_value_sq ,3) ) )
#     print("Done. pdf/molclock/" + subject_id + "_assembleBased.pdf")



def infer_non_snp_events(mutantAF,p,distance_for_nonsnp=500,fraction_covariance_nonsnp=0.98):
    ''' Mutations that  covary in close physical proximity (unexpected by chance) are likely non-independent mutations due to HGT event '''
    # return index of p with evidence for recombination (nonsnp) , boolean of length p with nonsnp == True
    p_nonsnps_events = np.array([])
    cov_m = np.cov(mutantAF)
    dist_p = p[1:]-p[:-1] # list of bp differences between elements of p (len(p)-1)
    proxy_pi = dist_p < distance_for_nonsnp
    proxy_pi = np.argwhere(proxy_pi) # index p and other p in close proximity
    for i,p_i in enumerate(proxy_pi):
        roi_min , roi_max = p[p_i] - distance_for_nonsnp , p[p_i] + distance_for_nonsnp # define region of interest with mutations, where we know at least 1 other mutation is in close proximity
        p_bool_roi = np.all( (p>roi_min , p<roi_max ),axis=0) # bool of p in proximity
        cov_p_i = cov_m[p_i,:][0,:] # cov of focal SNP p_i with all other p
        cov_p_bool = cov_p_i > (np.max(cov_p_i) * fraction_covariance_nonsnp) # bool of all SNPs that have x*max(cov) with SNP. Note: The focal SNP has ~highest cov with self
        if np.sum(cov_p_bool[p_bool_roi]) > 1: # if high cov of focal SNP with at least one more SNP in close proximity
            p_nonsnps_events = np.append( p_nonsnps_events, p[p_bool_roi][cov_p_bool[p_bool_roi]])
    p_nonsnps_events = np.unique(p_nonsnps_events)
    p_idx_nonsnp = np.array([np.where(p==i)[0][0] for i in p_nonsnps_events]) # get index of for nonsnps
    p_bool_nonsnp = np.zeros(p.shape) # build bool array of length p with nonsnps == True
    if p_idx_nonsnp.size > 0:
        p_bool_nonsnp[p_idx_nonsnp] = 1
    p_bool_nonsnp = p_bool_nonsnp.astype(bool)
    print(p_idx_nonsnp.size,'/',p.size,'elements of candpos evidence for recombination event.')
    return [p_idx_nonsnp,p_bool_nonsnp]





