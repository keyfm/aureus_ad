#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 2019

@author: fmk
"""

## To have all modules setup please run in conda environment: spyder4_full_env.yml

# =============================================================================
# Import modules
# =============================================================================
# %% apy module (contains all functions)

import sys,os,re,subprocess,gzip
import glob,pickle,subprocess,random,time,datetime,collections,warnings
import numpy as np
from scipy import stats
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from BCBio.GFF import GFFExaminer # pip install bcbio-gff
from BCBio import GFF # might be redundant
from Bio import Phylo
from Bio import AlignIO
from Bio.Data import CodonTable
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor #NJ tree

from collections import OrderedDict
from pylab import * #what is that

from matplotlib import rc
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.font_manager import FontProperties
import seaborn as sns

# Read analysis.py module:
SCRIPTS_DIRECTORY = "../modules" # path to analysispy_module.py
sys.path.insert(0, SCRIPTS_DIRECTORY)

import analysispy_module as apy

# %% Filter Parameter
# for finding fixed mutations between samples
filter_parameter_sample_across_sites = {\
                                        'min_average_coverage_to_include_sample': 8,
                                        'min_basecalls_to_include_sample': 0.1, # remove samples that have too many undefined base (ie. N). added this filter.
                                        'max_mean_MinorAF': 0.03, #0.04, NEW. Few isolates have high mean minorAF (6-8%), which causes excess wrong calls (misidentifeid secondary lineage in S8 and S14). )0.04 > 0;03 removes addiitonal 8 isolates; 0.04 cuts out outlier based on hist. see line 341. UPDATE. 0.03 is req for removal misplaced 067-CD5/066-N8
                                        }

filter_parameter_site_per_sample = {\
                                    'min_maf_for_call' : 0.85, #on individual samples, calls
                                    'min_cov_per_strand_for_call' : 2,  #on individual samples, calls
                                    'min_qual_for_call' : 30,  #on individual samples, calls
                                    }

filter_parameter_site_across_samples = {\
                                        'max_fraction_ambigious_samples' : 0.05, #0.25, #across samples per position. UPDATE: 0.05 to make more core genome like
                                        'min_median_coverage_position' : 3, #across samples per position
                                        }



# %% Labels
refgenome = 'Saureus_USA300FPR3757' # ref genome used in Key et al.

# %% Enviornment set up
refdir = '/To/Be/Defined/By/User/' # folder that contains ref genome folder. NOTE: Inside folder the fasta should be: 'genome.fasta'
ref_genome_folder = refdir + 'Reference_Genomes/' + refgenome # contains genome.fasta

# working directory
os.chdir('./across_patient_analysis')

## how far upstream of the nearest gene to annotate something a promoter
# mutation (not used if no annotation)
promotersize=250;


# %% Define outliers
outliers_all = array(['004-RD', '041-CD6', '027-N2', '027-N12', '027-N5', '027-N14','027-N16', '027-N18', '027-N19', '027-N11', '027-N15', '027-N17','027-N8', '027-N01', '027-N10', '027-N13', '027-N3', '027-N4','027-N6', '027-N7', '027-N9', '078-N8', '078-N9', '061-AD1','186-AI2', '023-AD1', '061-AD4', '061-N4', '071-CI10', '112-CI4','121-CI1'], dtype='<U8')

# %% Other variables
NTs = np.array(['A','T','C','G'],dtype=object) # NTs='ATCG'

# %% Load in specimen information
# parse specimen info to csv
specimen_log_file='../within_patient_analysis/aureus_specimen_log_v3.csv'; # contains for each kit# the pat# and visit#
file_specimen = open(specimen_log_file,"r")
SpecimenLog = {}
SpecimenLog["Patient"] = [];SpecimenLog["Visit"] = [];SpecimenLog["Kit"] = [];
for line in file_specimen:
    if not line.startswith("Patient"):
        line = line.strip().split(",")
        SpecimenLog["Patient"].append(line[0])
        SpecimenLog["Visit"].append(line[1])
        SpecimenLog["Kit"].append(line[2])


# =============================================================================
# BELOW NO FURTHER ADJUSTMENTS NECESSARY
# =============================================================================

# %% ANALYSIS

# %% 
# =============================================================================
# Load data from candidate_mutation_
# =============================================================================
# Import candidate_mutation_table
cmtFile = 'snakemake_raw_data_processing/case/2-candidate_mutation_table/candidate_mutation_table.pickle.gz' # generated via snakemake
[quals,p,counts,in_outgroup,sampleNames,indel_counter] = apy.read_candidate_mutation_table_pickle_gzip(cmtFile)


# %% continue processing data     
# % indel counter reduction to indel count
# indel_counter: The first statistic is the number of reads (at this position and in this sample) that 
# support an indel. The second statistics is the number of reads (at this position and in this sample) that support a deletion.
# for now we need only the count for indels overall 
indel_counter = indel_counter[:,0,:].transpose()
# indel counter >> 50% indel row:1; row2 >> deletion

# %% Save everything using more permanent sample names to avoid overwriting
sampleNames_all=np.asarray(sampleNames,dtype=object)
quals_all=-quals;
counts_all=counts;
coverage_all = counts_all.sum(axis=1).transpose() # axis=1 == rows; transpose needed > rows: pos and col: samples

# %% Assign patient & visit number to each sample; Parse specimen name to human readable
patients_all = np.empty(len(sampleNames_all), dtype=object) #np.zeros(shape=(1, len(sampleNames_all)), dtype=np.int) # row w/ 0's 
visits_all = np.zeros(len(sampleNames_all), dtype=np.int) # row w/ 0's
locations_all = np.zeros(len(sampleNames_all), dtype=np.int) # row w/ 0's
locations_abbreviations = [['RD'],['RDa'],['RI'],['CD','Cd'],['CI', 'Ci', 'Cl'],['CIa'],['N'],['AD'],['AI'],['B','C','NI']]; # there is a 'B5' and a 'C12', which is up to now obscure; NI in subj4 (041-NI)
locations_long_names = ['Right-Popliteal-Fossa','Right-Popliteal-Fossa-reseq','Left-Popliteal-Fossa','Right-Cubital-Fossa','Left-Cubital-Fossa','Left-Cubital-Fossa-reseq','Nare','Right-Forearm','Left-Forearm','unknown'];
    
for i in range(0,len(sampleNames_all)):
    kit_idx = SpecimenLog['Kit'].index( sampleNames_all[i][0:3] )
    patients_all[i] =  SpecimenLog['Patient'][kit_idx]
    visits_all[i] =  int(SpecimenLog['Visit'][kit_idx])
    currID = re.sub('-|[0-9]','',sampleNames[i]) # replace everything but location identifier
    for l in locations_abbreviations:
        if currID in l:
            locationmatch = locations_abbreviations.index(l) # get location of sublist that contains currID
    locations_all[i] = locationmatch

patients_sampled = np.unique(patients_all)
    
# %% Display samples the NOT fullfill min_average_coverage_to_include_sample: 
lowcovsamples = sampleNames_all[ np.any([coverage_all.mean(axis=0) < filter_parameter_sample_across_sites['min_average_coverage_to_include_sample'], np.isin(sampleNames_all, outliers_all,invert=False) ],axis=0) ]
print(lowcovsamples) # lowcov and outliers!
# 217/1787

# %%
# =============================================================================
#     Read in genome information
# =============================================================================
[chrStarts, genomeLength, scafNames] = apy.genomestats(ref_genome_folder);

# %%
# =============================================================================
#     Extract allele and frequencies
# =============================================================================
contig_positions = apy.p2chrpos(p,chrStarts) # 1col: chr, 2col: pos on chr; for all p
[maf, maNT, minorNT, minorAF] = apy.div_major_allele_freq(counts) 
# NOTE: function assumes first 8 rows in counts == 4nucl fwd&rev! watch out if extended counts used!
# NOTE: maf==0 -> no data;minorAF==0 -> no minor allele; NT number corresponds to index in NTs [ATCG]


# =============================================================================
#  Filter samples for coverage, outliers, minorAF
# =============================================================================
goodsamples =  np.all([coverage_all.mean(axis=0) >= filter_parameter_sample_across_sites['min_average_coverage_to_include_sample'], np.isin(sampleNames_all, outliers_all,invert=True), np.mean(minorAF,axis=0) < filter_parameter_sample_across_sites['max_mean_MinorAF'] ],axis=0) # , 


# =============================================================================
#  Adjust data matrices to goodsamples only
# =============================================================================
sampleNames = sampleNames_all[goodsamples]
counts = counts_all[goodsamples , : , : ] # keep only level (samples) that fullfil filter!
quals = quals_all[ : , goodsamples ]
maNT = maNT[:,goodsamples]
maf = maf[:,goodsamples]
minorNT = minorNT[:,goodsamples]
minorAF = minorAF[:,goodsamples]
coverage = coverage_all[ : ,goodsamples]
patients = patients_all[goodsamples]
visits = visits_all[goodsamples]
locations = locations_all[goodsamples]
indels = indel_counter[:,goodsamples]

num_samples = len(sampleNames)

coverage_forward_strand = counts[:,0:4,:].sum(axis=1).transpose()
coverage_reverse_strand = counts[:,4:8,:].sum(axis=1).transpose()
 

# %%
# =============================================================================
# Define outgroup nucleotide idx
# =============================================================================
refnt = apy.extract_outgroup_mutation_positions(ref_genome_folder, apy.p2chrpos(p,chrStarts));
print(np.unique(refnt)) # sanity check

# When no outgroup defined: refnt ~= ancnt:
ancnt = refnt

## build a pos x samples matrix for ancestral (ref) positions (based on NTs index)
## transform allele to nunmeric:0:3 for NTs->ATCG; all non-ATCG are "4"!; default value 9, should all be removed after loop!
ancnti_m = np.full(ancnt.shape, 9)
for idx,allele in enumerate(ancnt):
   if allele in NTs:
       ancnti_m[idx,] = np.where(NTs==allele)[0][0] # strip down to index number
   else:
       ancnti_m[idx,] = 4
ancnti_m = np.tile(ancnti_m,(num_samples,1)).transpose() # build 2D matrix



# %%
# =============================================================================
#     Find positions with fixed mutations
# =============================================================================
# Define mutations we do not trust in each and across samples.
# added indel filer! (less than 50% of reads had site are allowed to support indel)
calls = maNT
calls[ (quals < filter_parameter_site_per_sample['min_qual_for_call']) ] = 4
calls[ (maf < filter_parameter_site_per_sample['min_maf_for_call']) ] = 4
calls[ (coverage_forward_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call']) ] = 4
calls[ (coverage_reverse_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call']) ] = 4
calls[ (indels > (0.5*coverage) ) ] = 4 
# NOTE: reassigning calls to 4 changes also maNT!


siteFilt = np.any(( (calls>3).sum(axis=1) >= (num_samples * filter_parameter_site_across_samples['max_fraction_ambigious_samples']) \
                     ,np.median( coverage, axis=1) < filter_parameter_site_across_samples['min_median_coverage_position'] ),axis=0)
calls[ siteFilt ,:] = 4 # sites that fail qc -> 4, for all samples      
print('Done calls filter.')


# =============================================================================
#  Save mutQual, mutQualIsolates
# =============================================================================

# NOTE: func below takes forever with many SNPs...saved below
[mutQual, mutQualIsolates] = apy.ana_mutation_quality(calls,quals) # get FQ value for SNP across samples. mutQualIsolates contains sample indices for sample pair FQ based on. 
pickle_mutquals = [mutQual, mutQualIsolates]    

afile = open("mutQuals.py.pk1", 'wb')
pickle.dump(pickle_mutquals, afile)
afile.close()

# afile = open("mutQuals.py.pk1", 'rb')
# [mutQual, mutQualIsolates] = pickle.load(afile)
# afile.close()


## continue from above save
fixedmutation = (calls != ancnti_m) & (calls < 4) & (np.tile(mutQual,(1,num_samples)) >= 1) # nan in mutQual gives a warning >> fine (nan>False>excluded)
hasmutation = np.any( (fixedmutation == True , ), axis=0) # diversemutation == True #> not used
filteredpos = np.where( np.sum(hasmutation, axis=1) > 0 )[0] # NOTE: goodpos is INDEX of good positions for p!

goodpos = filteredpos #np.delete(filteredpos,non_variable_in_subject_idx)
print(str(len(goodpos)) + ' goodpos found.') 



# =============================================================================
#  Filter samples for minimum basecalls
# =============================================================================
## filter samples for min #bases called [used all 90k SNPS!]
filter_base_threshold_samples = len(goodpos)*(filter_parameter_sample_across_sites['min_basecalls_to_include_sample']) # 1-cutoff, because we test N not ACTG
samples_bool_w_minBasecalls = (calls != 4).sum(axis=0) >= filter_base_threshold_samples # get idx for samples with too little data


# %% Write sampleNames of goodSamples w/o outgroup to file
# Necessary to filter coverage matrix and remove false positives
# np.savetxt('sampleNames_goodsamples_noOut.txt',np.delete(sampleNames,outgroup_idx),fmt='%s')
np.save('sampleNames_goodsamples_noOut', sampleNames[samples_bool_w_minBasecalls]) # appends .npy


# =============================================================================
# Get multi-sample fasta for ML tree caclulation
# =============================================================================  

# define pos to use in tree
goodpos2useTree = goodpos #(1:1000); %trim for easier tree view; TDL called quality_positions

# define isolates to plot
samplestoplot = samples_bool_w_minBasecalls

# get data and filter pos
calls_for_treei = calls # equal to maNT which is propagated between both matrices   
calls_for_treei = calls_for_treei[ goodpos2useTree, : ] # numpy broadcasting of row_array and col_array requires np.ix_()

# reduce tree data to samplestoplot
calls_for_treei = calls_for_treei[ : , samplestoplot ] # numpy broadcasting of row_array and col_array requires np.ix_()

# translate index to nucleotide
calls_for_tree = apy.idx2nts(calls_for_treei) # ATCGN translation
# add reference nucleotide for all positions
outgroup_nts = apy.extract_outgroup_mutation_positions(ref_genome_folder, apy.p2chrpos(p[goodpos2useTree],chrStarts));

# build sampleNames w/ metainfo
treesampleNamesLong = apy.annotate_sampleNames(sampleNames[samplestoplot],locations_long_names,patients[samplestoplot],visits[samplestoplot],locations[samplestoplot]) 

# sampleNamesDnapars : max 10c. Use numeric with 10c (works for up to 10^10 samples! )
sampleNamesDnapars = [ "{:010d}".format(i) for i in range(len(samplestoplot))]

calls_for_tree = np.concatenate((outgroup_nts[:, None],calls_for_tree),axis=1) # first column now outgroup_nts; outgroup_nts[:, None] to make ndims (2) same for both
sampleNamesDnapars = np.append('Sref',sampleNamesDnapars) # add name for outgroup
treesampleNamesLong = np.append('Sref',treesampleNamesLong) # add name for outgroup

# build parsimony tree; added flag "buildTree": if false > only build dnaparse/fasta input files but do not calc tree 
apy.generate_tree(calls_for_tree,treesampleNamesLong,sampleNamesDnapars,refgenome,filetag='',buildTree=False)

#### > Generates a multi-fasta labeled YYYY-MM-DD-HH-MM-SS.fasta that can be used for phylogeny reconstruction.
taphAD/phylogen/2021_08_1700spl_USA300/run2_partialDel005_minorAF003 .



# =============================================================================
# Pairwise differences
# =============================================================================  
# Assess quantitative differnce between lineages and double colonizers >> Get pairwise differences for all isolates
# USA300 alignment > all goodpos where both isolates have ACTG.
# Build martix with all pairs
# Pickle matrix, num p usable sites, as well as sample Names in dict pairwise_data


## > Use calls_for_treei in dbl loop to get pairwise differences. does not contain SREF
pairwise_diff = np.array([],dtype=int)

for i in arange(samplestoplot.shape[0]):
    print(i)
    for j in arange(samplestoplot.shape[0]):
        if i == j: # i j same > no diff
            pairwise_diff = np.append(pairwise_diff, 0)
        elif j < i: # j<i : i vs. j alerady calculated
            idx_first_calculation = j*samplestoplot.shape[0]+i
            pairwise_diff = np.append(pairwise_diff, pairwise_diff[idx_first_calculation])
        else:
            bool_diff = np.not_equal(calls_for_treei[:,i], calls_for_treei[:,j]) # get diffs (True) for all pairwise diffs
            bool_i_j_nonN = np.logical_and( calls_for_treei[:,i] != 4  , calls_for_treei[:,j] != 4) # get bool (True) of pos that are ACTG and not N in BOTH isolates (i and j)
            pairwise_diff = np.append(pairwise_diff, np.sum( bool_diff[bool_i_j_nonN] ) ) # count all True which equals num diffs
            
pairwise_diff = pairwise_diff.reshape( (len(samplestoplot), len(samplestoplot)) )


# get list with # of usable (non-4) calls
p_usable = np.array([],dtype=int)
for i in arange(samplestoplot.shape[0]):
    bool_p_usable = calls_for_treei[:,i] != 4
    p_usable = np.append(p_usable, np.sum( bool_p_usable ) )


## pickle

pairwise_data = {}
pairwise_data['treesampleNamesLong'] =  treesampleNamesLong[1:] # remove SRef
pairwise_data['pairwise_diff'] =  pairwise_diff
pairwise_data['p_usable'] =  p_usable


with open('pairwise_data.pickle', 'wb') as handle:
    pickle.dump(pairwise_data, handle, protocol=pickle.HIGHEST_PROTOCOL)



 
    