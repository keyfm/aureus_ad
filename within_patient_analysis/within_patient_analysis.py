#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 2020

@author: fmk
"""

## To have all modules setup please run in conda environment: spyder4_full_env.yml

# %%
# =============================================================================
# SNP-based analyses:
# =============================================================================
# - parsimony tree
# - SNP-specific tree coloring
# - barplots fwd/rev for SNPs
# - parallel evolution analysis
# - molecular clock analysis
# - read in unified ortholog table (creation NOT yet part of assembly snakemake)
# - allele frequency change analysis
# - read-in fucntions from apy (analysis.py) module

# %% apy module (contains all functions)

import sys,os,re,subprocess,gzip,json
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
from pylab import * 

from matplotlib import rc
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.font_manager import FontProperties
import seaborn as sns


SCRIPTS_DIRECTORY = "../modules"
sys.path.insert(0, SCRIPTS_DIRECTORY)

import analysispy_module as apy

# import importlib
# importlib.reload(apy)

# %% Filter Parameter -- IMPORTANT TO ADJUST !!!!
# Much of this will need to vary according to your reference genome,
# coverage, and particular samples
# consider to add SD filter in addition to coverage. Isolates of other species that harbor plasmid present in pangenome at very high cov can otherwise meet the coverage threshold even though FP

# for finding fixed mutations between samples
filter_parameter_sample_across_sites = {\
                                        'min_average_coverage_to_include_sample': 8, 
                                        'min_basecalls_to_include_sample': 0.1, # remove samples that have too many undefined base (ie. N). added this filter.
                                        }

filter_parameter_site_per_sample = {\
                                    'min_maf_for_call' : 0.85, #on individual samples, calls
                                    'min_cov_per_strand_for_call' : 2,  # on individual samples, calls
                                    'min_qual_for_call' : 30,  #on individual samples, calls
                                    }

filter_parameter_site_across_samples = {\
                                        'max_fraction_ambigious_samples' : 0.25, #across samples per position
                                        'min_median_coverage_position' : 3, #across samples per position
                                        'fraction_covariance_of_snps_across_isolates_for_removal' :0.90, # covariance of SNPs in close proximity (500b) signals recombination and will be removed
                                        }


## how far upstream of the nearest gene to annotate something a promoter
## mutation (not used if no annotation)
promotersize=250;

# %% Other variables
NTs = np.array(['A','T','C','G'],dtype=object) # NTs='ATCG'

# %% Load in specimen information
# parse specimen info to csv
specimen_log_file='aureus_specimen_log_v3.csv'; # updated for new samples
file_specimen = open(specimen_log_file,"r")
SpecimenLog = {}
SpecimenLog["Patient"] = [];SpecimenLog["Visit"] = [];SpecimenLog["Kit"] = [];
for line in file_specimen:
    if not line.startswith("Patient"):
        line = line.strip().split(",")
        SpecimenLog["Patient"].append(line[0])
        SpecimenLog["Visit"].append(line[1])
        SpecimenLog["Kit"].append(line[2])

# %% Read Date and Timing info per Patient
pat_visit_timing = pd.read_csv('../metadata/patient_visit_date_timepast.csv')

# %% Read in gene ortholog information; build during assembly snakemake step
ortholog_data = 'snakemake_raw_data_processing/assembly/6-ortholog_identification/annotation_orthologs.tsv' 
ortholog_df = pd.read_csv(ortholog_data,sep="\t",index_col=0)
        
# %% Read subject-specific outgroups
# isolate with highest coverage from closest lineage (as in Figure 1)
outgroup_specifier = pd.read_csv('outgroups_subject_top_cov.csv') 

             
# %% Get all subject folder names (used in working directory and ref genome folder)
subject_fld_name = os.popen('ls snakemake_raw_data_processing/case/').read().split('\n')
subject_fld_name = np.array([i for i in subject_fld_name if i != ""]) # turn to np.array and remove empty value


# %% Define outliers for each subject
## minor lineages and cross contaminants
outliers_all = np.array(['067-AD5','069-RD9','069-RI1','069-RI2','069-RI4','069-RI7','121-CD1','121-CD10','121-CD2','121-CD3','121-CD4','121-CD5','121-CD6','121-CD7','121-CD8','121-CD9','121-CI10','121-CI3','121-CI8','121-CI9','121-RD1','121-RD3','121-RD4','121-RD5','121-RD6','121-RD7','121-RD8','121-RD9','121-RI10','121-RI2','121-RI3','121-RI4','121-RI5','121-RI6','121-RI9'])


# %% TIMESTAMP
# used in output files
ts = time.time()
timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S') # WARN: fixed timestamp not yet build in
 

# %% Define variables to store subject-specific results 
para_evo_cand = np.array([],dtype=object)
mut_genes_df = {}
snp_freq_shifts_anno_lod = []
regression_res_lod = [] # stores all the results of the regression, and turned to pandas DF and saved
annotation_mutation_allParaSignal = {} # store all annotation_mutation for all detected candidates
mol_clock_data_dc = {} # dc to store for each subject the calculated data for mol clock estimation
allP_gene_paraevo = {} # store all para evo candidates that fullfill bonferroni

# %% Initiate loop over all subjects
for subj_idx,subject_fld_label in enumerate(subject_fld_name):
    print('\n'+'Process ' + subject_fld_label )

    # %% Labels
    refgenome = subject_fld_label

    # cd(workingdir)
    os.mkdir('analysis/'+refgenome)
    os.chdir('analysis/'+refgenome)
        
    # %% Enviornment set up -- probably won't need to change# 
    # requires a 
    ref_genome_folder = 'snakemake_raw_data_processing/assembly/3-spades/' + refgenome 
    

    # %% 
    # =============================================================================
    # Load data from candidate_mutation_
    
    # =============================================================================
    # load('candidate_mutation_table')
    # 
    # Import candidate_mutation_table
    cmtFile = '../../snakemake_raw_data_processing/case/2-candidate_mutation_table/candidate_mutation_table.pickle.gz'
    [quals,p,counts,in_outgroup,sampleNames,indel_counter] = apy.read_candidate_mutation_table_pickle_gzip(cmtFile)

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

    # %%
    # =============================================================================
    #     Read in genome information
    # =============================================================================
    [chrStarts, genomeLength, scafNames] = apy.genomestats(ref_genome_folder);

    
    # %% Display samples the NOT fullfill min_average_coverage_to_include_sample: 
    lowcovsamples = sampleNames_all[ np.mean(coverage_all, axis=0) < filter_parameter_sample_across_sites['min_average_coverage_to_include_sample'] ]
    print(lowcovsamples)
    if len([i for i in lowcovsamples if re.search("_o",i)]) > 4 or subject_fld_label == 'StaphAureus_fmk_2-control':
        print('Outgroup filtered due to low coverage >> skip')
        continue
    
   
    # %% 
    # =============================================================================
    #     Remove undesired samples based on outlier dict &| coverage &| genome-wide coverage
    # =============================================================================
    # read in cov file. #num genome-wide bases cov>=8, generated during SM
    numBasesCovGenome = pd.read_csv('numPos_covThreshold8.csv') 
    minNumBases = genomeLength * filter_parameter_sample_across_sites['min_basecalls_to_include_sample']

    # %% Define goodsamples and filter data
    goodsamples =  np.all([coverage_all.mean(axis=0) >= filter_parameter_sample_across_sites['min_average_coverage_to_include_sample'], np.isin(sampleNames_all, outliers_all,invert=True) , np.array(numBasesCovGenome['numpos_covthreshold']) > minNumBases],axis=0) # , 
            
    sampleNames = sampleNames_all[goodsamples]
    counts = counts_all[goodsamples , : , : ] # keep only level (samples) that fullfil filter!
    quals = quals_all[ : , goodsamples ]
    coverage = coverage_all[ : ,goodsamples]
    patients = patients_all[goodsamples]
    visits = visits_all[goodsamples]
    locations = locations_all[goodsamples]
    indels = indel_counter[:,goodsamples]
    
    num_samples = len(sampleNames)
    
    coverage_forward_strand = counts[:,0:4,:].sum(axis=1).transpose()
    coverage_reverse_strand = counts[:,4:8,:].sum(axis=1).transpose()
    
    ## record all goodsamples
    np.savetxt('goodsamples.txt',apy.annotate_sampleNames(sampleNames,locations_long_names,patients,visits,locations),fmt="%s")
   
    # %% Breakpoint: Too few samples passed filter
    if np.sum(goodsamples) < 2:
        print("Too few samples fullfill filter criteria! >> skip: " + refgenome)
        continue

    # %%
    # =============================================================================
    # Extract refnt and define out/in-group bools
    # =============================================================================
    refnt = apy.extract_outgroup_mutation_positions(ref_genome_folder, apy.p2chrpos(p,chrStarts));
    refnti = apy.nts2idx(refnt)
    refnti_m = np.tile(refnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele
    # print(np.unique(refnt)) # sanity check
    
    # When no outgroup defined: refnt ~= ancnt:
    #ancnt = refnt   
    
    ## Estimate outgroup (ancestral allele) from ALL samples added as outgroup to SM pipeline (ancnti* == major allele)
    pattern = re.compile('_o') # string tag in sample Name to identify outgroup in staphAD data
    outgroup_name = np.array(list(filter(pattern.search, list(sampleNames))))
    outgroup_bool = np.in1d(sampleNames , outgroup_name)
    
    # ingroup array (bool, idx) used later
    ingroup_bool = np.invert(outgroup_bool)
    ingroup_idx = np.nonzero(ingroup_bool)[0]

    ## bool for all single outgroup sample AND ingroup-samples with dedicated single outgroup sample
    # based Fig1-USA300 tree we select the phylogentically closest (w/ highest cov vs.USA300) sample as single outgroup for tree
    # NOTE: patients w/ _a /_b lineages will be assigned the outgroup sample from _a. That is correct for the 5 important subjects and should also be correct for all other as we do not plot the _b lineage
    outgroup_currSub = outgroup_specifier['outgroup_sample'][ outgroup_specifier['patient'] == subject_fld_label.split('_')[-1] ].to_numpy()
    if outgroup_currSub.size == 1:
        # single lineage for subject
        outgroup_spl_idx = [i for i,item in enumerate(sampleNames) if outgroup_currSub[0] in item]
        ingroup_wOutSpl_bool = np.copy(ingroup_bool)
        ingroup_wOutSpl_bool[outgroup_spl_idx] = True
    elif outgroup_currSub.size == 2:
        # a/b lineage for subject
        idx = [i for i,item in enumerate(outgroup_specifier['patient_clade'][ outgroup_specifier['patient'] == subject_fld_label.split('_')[-1] ].to_numpy()) if item == subject_fld_label.split('_')[-1]+"_a"] # index in outgroup_currSub that corresponds to _a lineage
        outgroup_spl_idx = [i for i,item in enumerate(sampleNames) if outgroup_currSub[idx][0] in item]
        ingroup_wOutSpl_bool = np.copy(ingroup_bool)
        ingroup_wOutSpl_bool[outgroup_spl_idx] = True
    else:
        # outgroup_currSub should only be of size 1 or 2 otherwise something is wrong and requires adjustment.
        warnings.warn("Warning: outgroup sample remains undefined and will not be included in the tree!")
        ingroup_wOutSpl_bool = np.copy(ingroup_bool)
    
    ingroup_wOutSpl_onlyIngroup_bool = np.array([False if "_o" in s else True for s in sampleNames[ingroup_wOutSpl_bool]])
   
    
    # %%
    # =============================================================================
    #     Extract allele and frequencies
    # =============================================================================
    contig_positions = apy.p2chrpos(p,chrStarts) # 1col: chr, 2col: pos on chr; for all p
    [maf, maNT, minorNT, minorAF] = apy.div_major_allele_freq(counts) 
    # NOTE: function assumes first 8 rows in counts == 4nucl fwd&rev! watch out if extended counts used!
    # NOTE: maf==0 -> no data;minorAF==0 -> no minor allele/ or no major allele; NT number corresponds to index in NTs [ATCG] or if maf==0 > NA == 4  

    
    # %%
    # =============================================================================
    # Make some basic structures for finding mutations
    # =============================================================================
    mutantAF = np.zeros(maNT.shape)
    mutantAF[maNT != refnti_m] = maf[ maNT != refnti_m];        

    
    # %%
    # =============================================================================
    #     Find positions with fixed mutations
    # =============================================================================
    # Define mutations we do not trust in each and across samples.
    # goodpos are indices of p that we trust
    
    ## Filter per mutation
    # added indel filer! (less than 50% of reads had site are allowed to support indel)
    calls = maNT
    calls[ (quals < filter_parameter_site_per_sample['min_qual_for_call']) ] = 4
    calls[ (maf < filter_parameter_site_per_sample['min_maf_for_call']) ] = 4
    calls[ (coverage_forward_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call']) ] = 4
    calls[ (coverage_reverse_strand < filter_parameter_site_per_sample['min_cov_per_strand_for_call']) ] = 4
    calls[ (indels > (0.5*coverage) ) ] = 4 
    
    ## Filter per site across samples
    # Ignore here outgroup samples!
    siteFilt = np.any(( (calls[:,ingroup_bool]>3).sum(axis=1) >= ((num_samples-np.sum(outgroup_bool)) * filter_parameter_site_across_samples['max_fraction_ambigious_samples']) \
                         ,np.median( coverage[:,ingroup_bool], axis=1) < filter_parameter_site_across_samples['min_median_coverage_position'] ),axis=0)
    calls[ siteFilt ,:] = 4 # sites that fail qc -> 4, for all samples incl. outgroup     
        
    # NOTE: func below takes forever with many SNPs...saved below
    [mutQual, mutQualIsolates] = apy.ana_mutation_quality(calls[:,ingroup_bool],quals[:,ingroup_bool]) # get FQ value for SNP across samples. mutQualIsolates contains sample indices for sample pair FQ based on. 
    mutQual = np.nan_to_num(mutQual, nan=-1) # turn mutQual nan's to -1; necessary to avoid later warning
    
    ## translate filtered calls of ingroup into goodpos. mutations we believe. fixedmutation part removed in v6.
    hasmutation = (calls != refnti_m) & (calls < 4) & (np.tile(mutQual,(1,num_samples)) >= 1) # consider only ingroup samples; mutQual >= 1 is very loose. Important filter with low qual data! refnt not ancnt!!!
    hasmutation[:,outgroup_bool] = False # put outgroup samples 4 in order to identify ingroup mutations only
    
    candpos = np.where( np.sum(hasmutation, axis=1) > 0 )[0] # NOTE: candpos/goodpos is INDEX of good positions for p!

    ## identify mutations due to recombination
    # defined as mutation that has its max(covarying) mutation in close physical proximity
    # Based on ingroup only. Important as outgroup samples can interfere when covarying together with a single nonsnp but not the entire recombinatory region  
    # added to mutation filter below at hasmutation build up
    # do candpos first to reduce number of SNPs to analyse...computational cost high if >40k SNPs
    [nonsnp_p_idx,nonsnp_p_bool] = apy.infer_non_snp_events( mutantAF[ np.ix_(candpos, ingroup_bool) ],p[candpos],distance_for_nonsnp=500,fraction_covariance_nonsnp=filter_parameter_site_across_samples['fraction_covariance_of_snps_across_isolates_for_removal']) # original setting: 0.98 # np.ix_ required to slice 2D array with index (rows) and boolean ()
    
    ## remove mutations within a gene if only present in single isolate    
    goodpos = candpos[ np.invert(nonsnp_p_bool) ] # NOTE: candpos/goodpos is INDEX of good positions for p! 

    print(goodpos.size,'goodpos found.')
    
   
    # %% 
    # =============================================================================
    #     Breakpoint: Too few positions pass filter
    # =============================================================================
    if len(goodpos) < 2:
        print("Too few positions after filter! >> skip: " + refgenome)
        continue

    calls_outgroup = calls[:,outgroup_bool] 
    ancnti = apy.major_allele(calls_outgroup) # NOTE: the filter criteria (cov,qual etc.) are not applied before major allele call
    ancnti_m = np.tile(ancnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele    

    ## outsplnti == phylogenetically closest high-cov outgroup
    if subject_fld_label == 'subject_26':
        ## reroot tree for pat26 dominant lineage according to treewas [see Methods]
        calls[goodpos[92],outgroup_spl_idx]=0 # W146C
        calls[goodpos[96],outgroup_spl_idx]=1 # A104T
        calls[goodpos[119],outgroup_spl_idx]=3 # K108R
        calls[goodpos[160],outgroup_spl_idx]=0 # V21I
        
    outsplnti = calls[:,outgroup_spl_idx][:,0]
    outsplnti_m = np.tile(outsplnti,(num_samples,1)).transpose() # build 2D matrix with outgroup (ancestral) allele    
    print('Unresolved goodpos in outsplnti outgroup:',np.sum(outsplnti[goodpos]==4))

    ## out chimera. based on outsplnti but all NA (==4) are replaced by major allele call in ancnti
    # we use this for molecular clock and derived allele freq (but not tree!)
    outchimerasplancnti = np.array([ancnti[i] if nti==4 else nti for i,nti in enumerate(outsplnti) ])
    print('Unresolved goodpos in chimera outgroup:',np.sum(outchimerasplancnti[goodpos]==4))
    
    # %%
    # =============================================================================
    #  Display table (with or without annotation)
    # =============================================================================
    goodpos2use = goodpos; # leave. downsize only goodpos2useTree below if necessary
    order = np.arange(sampleNames.shape[0])
    num_contigs = np.max(contig_positions[:,0]);
    contig_lengths = np.append(chrStarts[1:], genomeLength) - chrStarts
    
    # Uncomment the following if there is an annotated genome
    # NOTE: annotation is read from *.gff file in ref_folder! Prokka builds gff. Otherwise NCBI provides gff for available genomes.
    annotation_genes = apy.parse_gff(ref_genome_folder,scafNames,ortholog_df[subject_fld_label],forceReDo=False) # ref_genome_folder+"/annotation_genes.pandas.py.pk1"
    annotation_mutations = apy.annotate_mutations(annotation_genes , p[goodpos2use] , refnti_m[goodpos2use,:] , outsplnti_m[goodpos2use,:] , calls[goodpos2use,:] , counts[:,:,goodpos2use] , hasmutation[goodpos2use,:], mutQual[goodpos2use,].flatten() , promotersize , ref_genome_folder) # extract relevant annotation info for each SNP    
    
    # %% 
    # =============================================================================
    #     parsimony tree. (incl. samples filter)
    # =============================================================================    
    # define pos to use in tree
    goodpos2useTree = goodpos #(1:1000); %trim is > 1000 SNPs (too computational intensive for desktop/laptop)
    
    # get data and filter for goodpos
    calls_for_treei = calls; 
    calls_for_treei = calls_for_treei[ goodpos2useTree, : ]
    calls_for_treei = calls_for_treei[ : , ingroup_wOutSpl_bool ]
    
    # build sampleNames w/ metainfo
    treesampleNamesLong = apy.annotate_sampleNames(sampleNames,locations_long_names,patients,visits,locations) # all sample names translated
    treesampleNamesLong = treesampleNamesLong[ingroup_wOutSpl_bool] # remove outgroup samples
        
    # sampleNamesDnapars : max 10c. Use numeric with 10c (works for up to 10^10 samples! )
    sampleNamesDnapars = [ "{:010d}".format(i) for i in range(len(treesampleNamesLong))]

    # translate index to nucleotide
    calls_for_tree = apy.idx2nts(calls_for_treei) # ATCGN translation
    # add reference nucleotide for all positions
    refgenome_nts = apy.extract_outgroup_mutation_positions(ref_genome_folder, apy.p2chrpos(p[goodpos2useTree],chrStarts));
    calls_for_tree = np.concatenate((refgenome_nts[:, None],calls_for_tree),axis=1) # first column now refgenome_nts; refgenome_nts[:, None] to make ndims (2) same for both

    sampleNamesDnapars = np.append(['Sref'],sampleNamesDnapars) # add name for outgroup
    treesampleNamesLong = np.append(['Sref'],treesampleNamesLong) # add name for outgroup

    # ## build parsimony tree; added flag "buildTree": if false > only build dnaparse/fasta input files but do not calc tree 
    apy.generate_tree(calls_for_tree,treesampleNamesLong,sampleNamesDnapars,refgenome,filetag=refgenome,buildTree='PS')
    
    # %% Write sampleNames of goodSamples w/o outgroup to file
    # Necessary to filter coverage matrix and remove false positives
    np.savetxt('sampleNames_goodsamples_noOut.txt',sampleNames[ingroup_bool],fmt='%s')
    np.save('sampleNames_goodsamples_noOut', sampleNames[ingroup_bool]) # appends .npy
    np.savetxt('sampleNames_goodsamples_treeNamesLongAll.txt',treesampleNamesLong,fmt='%s')
    
    # %%
    # =============================================================================
    #     Make a tree for each SNP location 
    # =============================================================================
    
    ## create empty tree_counting folder and add for_tree_labeling.csv
    apy.build_table_for_tree_labeling(apy.p2chrpos(p[goodpos2use], chrStarts),treesampleNamesLong,calls_for_tree)
    
    ## add tree for each mutation dsiplaying mutation (countMutations.py)
    os.chdir('tree_counting')
    apy.countMutations("../"+refgenome+"_latest.nwk.tree","for_tree_labeling.csv")
    os.chdir('../')
    
    # =============================================================================
    # %% Store SNP table
    # =============================================================================
    ## SOM SNP Table. Should be a function.
    ## get NT,sampleName df
    # build sampleNames w/ metainfo
    sampleNamesLong = apy.annotate_sampleNames(sampleNames,locations_long_names,patients,visits,locations) # all sample names translated
    sampleNamesLong = sampleNamesLong[ingroup_wOutSpl_bool] # remove outgroup samples
    sampleNamesLong = sampleNamesLong[ingroup_wOutSpl_onlyIngroup_bool] # remove outgroup
    
    # translate index to nucleotide
    calls_for_treei = calls; 
    calls_for_treei = calls_for_treei[ goodpos2useTree, : ]
    calls_for_treei = calls_for_treei[ : , ingroup_wOutSpl_bool ]
    calls_for_treei_ingroup = calls_for_treei[ : , ingroup_wOutSpl_onlyIngroup_bool ]
    calls_for_tree = apy.idx2nts(calls_for_treei_ingroup) # ATCGN translation
    snp_data = pd.DataFrame(calls_for_tree,columns=sampleNamesLong)

    ## get snp metadata
    # Contig	Location	Cluster ID	NCTC9343 homolog	Assembly locus	Gene	Protein annotation (Prokka)	Type*	Ancestor allele	Amino acid changes	Nucleotide position in the gene**
    snp_metadata = annotation_mutations.copy()
    # for I/P mutations: turn all to I (P somewhat random); add up/downstream info to 0tag and gene description to product
    for i,row in snp_metadata.iterrows(): # I/P SNV
        if row.isnull()['locustag']:
            if row.isnull()['locustag1'] and not row.isnull()['locustag2']: # no preceding gene. start of contig
                snp_metadata.at[i,'locustag'] = "NA;"+str(int(row['distance2']))+":"+row['locustag2']
                snp_metadata.at[i,'gene'] = "NA;"+row['gene2']
                snp_metadata.at[i,'product'] = "NA;"+row['product2']
                snp_metadata.at[i,'type'] = 'I'
            elif row.isnull()['locustag2'] and not row.isnull()['locustag1']: # no subsequent gene. end of contig
                snp_metadata.at[i,'locustag'] = str(row['distance1'])+":"+str(row['locustag1'])+";NA"
                snp_metadata.at[i,'gene'] = row['gene1']+";NA"
                snp_metadata.at[i,'product'] = row['product1']+";NA"
                snp_metadata.at[i,'type'] = 'I'
            elif row.isnull()['locustag1'] and row.isnull()['locustag2']: # no annotated gene on contig
                snp_metadata.at[i,'locustag'] = "NA"
                snp_metadata.at[i,'gene'] = "NA"
                snp_metadata.at[i,'product'] = "NA"
                snp_metadata.at[i,'type'] = 'I'
            else: # intergenic SNV with preceding/subsequent gene                
                snp_metadata.at[i,'locustag'] = str(row['distance1'])+":"+str(row['locustag1'])+";"+str(int(row['distance2']))+":"+row['locustag2']
                snp_metadata.at[i,'gene'] = row['gene1']+";"+row['gene2']
                snp_metadata.at[i,'product'] = row['product1']+";"+row['product2']
                snp_metadata.at[i,'type'] = 'I'

    snp_metadata = snp_metadata[['chr','pos','type','muts','locustag','gene','loc1','loc2','strand','product','nt_pos','aa_pos','nt_ref','nt_alt','nt_anc']]
    snp_metadata[['nt_anc']] = snp_metadata[['nt_anc']].replace('.','?') # turn NA to '?', similar to NT data 
    snp_table = pd.concat([snp_metadata,snp_data.reset_index(drop=True)], axis=1)

    ## store
    with open('snp_table.csv', 'w') as file:
        snp_table.to_csv(file, header=True, index=False)


    # %%
    # =============================================================================
    #     Plot barchart with fwd/rev read count for each position for each allele
    # =============================================================================
    # get dataframes that carry fwd/rev coverage for all bases (incl. N)      
    [lod_fwd_cov,lod_rev_cov] = apy.build_dataframe_coverage_info(goodpos2useTree,NTs,sampleNames,maNT,minorNT,coverage_forward_strand,coverage_reverse_strand,maf,minorAF)
    
    # get chr and pos for reelvant SNPs
    chr_pos_gp = contig_positions[goodpos2useTree,]
    
    # loop over SNPs and plot. 
    # All results in: pdf/coverage_snp_fwd_rev/chr_poos_locID_anno_qual.pdf
    apy.plot_coverage_fwd_rev_stacked(chr_pos_gp,annotation_mutations,lod_fwd_cov,lod_rev_cov,timestamp)
    
    
    # %%
    # =============================================================================
    # Estimate substitution rate
    # =============================================================================
    ## DONE using all subjects mapped to USA300
    allmuts = np.array(['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG'],dtype=object)
    # AT, TA  0
    # AC, TG  1
    # AG, TC  2
    # GC, CG  3
    # GT, CA  4
    # GA, CT  5
    allmuts_types = np.array([0,2,1,0,1,2,5,4,3,4,5,3])
    mutationalspectrum = [1/12] * 12 # uniform distribution. 
    
    ## loop over all ref/anc nucleotides and count mutation types. count twice for each direction
    allmuts = np.array(['AT','AG','AC','TA','TG','TC','GA','GT','GC','CA','CT','CG'],dtype=object)
    allmuts_counts = np.zeros( allmuts.shape,dtype=int )
    for i,ref in enumerate(calls_for_tree[:,0]): # Note: 1st col is reference!
        obs_allele = np.unique(calls_for_tree[i,1:])
        obs_allele = obs_allele[ ~(obs_allele == '?') & ~(obs_allele == ref) ]
        for j in obs_allele:
            idx = np.where( (ref+j)==allmuts)
            allmuts_counts[idx] += 1
            idx = np.where( (j+ref)==allmuts)
            allmuts_counts[idx] += 1
    mutationalspectrum = allmuts_counts/len(calls_for_tree[:,0])/2 # true mutationalspectrum
    # store
    afile = open("mutationalspectrum.py.pk1", 'wb')
    pickle.dump(mutationalspectrum, afile)
    afile.close()
    
    # %%
    # =============================================================================
    # Read in number of positions covered across genome
    # =============================================================================
    # obtain isolate specific-count of covered based required to correct observed mutation count (ie. for molecular clock inference and paraevo-poisson)
    # NOTE: This data has to be generated by user as it dependes on the assembly build. Here I put a placeholder
    numBaseGenome_covThreshold = 2800000

    # %%
    # =============================================================================
    # Parallel evolution module
    # =============================================================================
    # define dict with parameters for inference
    parameters = {
            'NumTrialsSim':1000,
            'Min_num_mutations_cand':2, # minimum number of mutations per gene candidate
            'Min_mutation_density_cand':1/1000, # minimum number of mutations per 1000 bp per gene candidate
            'ref_genome_folder':ref_genome_folder,
            'subjectID': refgenome, # used for naming pdf in pdf/adaptive_evo
            'substitution_spectrum':None, # put None if not calculated. Alternatively path
            'max_muts_per_gene_to_track': 100,
            'timestamp':timestamp
            }
    
    [res_cand_nummut,annotation_mutation_paraSignal] = apy.parallel_evo_module( goodpos2use, contig_positions , annotation_mutations , annotation_genes , parameters )

    # %%
    # =============================================================================
    # Report para evo results w/ Bonferroni Poisson P val 
    # =============================================================================
    # get poisson p value for observed candidates of para evo 

    ## variables for lambda (expectation) calculation
    tot_num_mutations = len(goodpos)
    median_genome_length = np.median(numBaseGenome_covThreshold['numpos_covthreshold'])

    ## bonferroni corrected p for subject (based on total number of annotated genes)
    num_annotated_genes = sum([t.shape[0] for t in annotation_genes if not t.empty])
    alpha_sig_level = 0.05 # pval considered significant for multiple testing correction
    p_bonf_corrected = alpha_sig_level/num_annotated_genes

    
    lod_gene_paraevo = []    
    for i in range(res_cand_nummut.shape[0]):
        para_gene_dict = {}
        ### calculate poisson lambda - expectation
        ## get gene of para evo candidate gene
        idxGeneAll = np.where(annotation_mutations['locustag']==res_cand_nummut[i,1])[0] 
        idxGene = idxGeneAll[0] # only take idx of first instance of gene
        lengthGene = len(annotation_mutations['sequence'][idxGene])
        my_lambda = tot_num_mutations * (lengthGene/median_genome_length)
        ### poisson cdf for num-obs_mut or more
        obs_mut = res_cand_nummut[i,2]
        p_poisson_obs = 1-stats.poisson.cdf(obs_mut-1,my_lambda)
        ### record if below bonferroni corrected p for that patient: p/#genes
        #if p_poisson_obs < p_bonf_corrected:
        # record all candidates
        para_gene_dict['Patient'] = subject_fld_label
        para_gene_dict['Contig'] = annotation_mutations['chr'][idxGene]
        para_gene_dict['Gene_Start'] = annotation_mutations['loc1'][idxGene]
        para_gene_dict['Gene_End'] = annotation_mutations['loc2'][idxGene]
        para_gene_dict['Gene_length'] = annotation_mutations['loc2'][idxGene]- annotation_mutations['loc1'][idxGene]
        para_gene_dict['Strand'] = annotation_mutations['strand'][idxGene]
        para_gene_dict['Gene_id_prokka'] = annotation_mutations['locustag'][idxGene]
        para_gene_dict['Gene_id'] = annotation_mutations['gene'][idxGene]
        para_gene_dict['Orthologtag'] = annotation_mutations['orthologtag'][idxGene]
        para_gene_dict['Product_description_prokka'] = annotation_mutations['product'][idxGene]
        para_gene_dict['Gene_id_orthologtag'] = annotation_mutations['orthologtag'][idxGene]
        para_gene_dict['Poisson_p_value'] = p_poisson_obs
        para_gene_dict['Bonferroni_p_value'] = p_bonf_corrected
        para_gene_dict['Numer_mutations'] = res_cand_nummut[i,2]
        para_gene_dict['ProbSim'] = res_cand_nummut[i,3]
        para_gene_dict['Mutations'] = [annotation_mutations['muts'][m][0] for m in idxGeneAll]
        para_gene_dict['Mutation_type'] = [annotation_mutations['type'][m] for m in idxGeneAll]
        para_gene_dict['Mutation_position'] = [annotation_mutations['pos'][m] for m in idxGeneAll]
        para_gene_dict['Translation'] = annotation_mutations['translation'][idxGene]
        para_gene_dict['Sequence'] = annotation_mutations['sequence'][idxGene]            
        lod_gene_paraevo.append(para_gene_dict)
    df_gene_paraevo = pd.DataFrame(lod_gene_paraevo)  # turn patient results to dataframe...
    if not df_gene_paraevo.empty:
        allP_gene_paraevo[subject_fld_label] = df_gene_paraevo # ... and store

           

    # %%    
    # =============================================================================
    #     Store candidate info and annotation mutations
    # =============================================================================
    # store only cand from subject with clonal strain colonisation, however, keep anno mutation for all

    if len(res_cand_nummut) > 0 and subject_fld_label not in ['subject_10','subject_21']:
        para_evo_cand = np.append(para_evo_cand,res_cand_nummut[:,0])
    mut_genes_df[subject_fld_label] = annotation_mutations # store all mutations
    
    if not annotation_mutation_paraSignal.empty:
        annotation_mutation_allParaSignal[subject_fld_label] = annotation_mutation_paraSignal
        
    # %%
    # =============================================================================
    #  Molecular Clock Module
    # =============================================================================
    # mutations polarized using outchimerasplanc which is based on phylogenetically close outgroup sample, with possible unknown bases (==4) filled up by major allele across all outgroups
    treesampleNamesLongNoRef = np.delete( treesampleNamesLong , np.where(treesampleNamesLong=='Sref')[0][0] ) # remove 'Sref' as it is not part of ingroup_bool
    list_visit_idx = apy.get_index_visit_sampleNames(visits[ingroup_bool],treesampleNamesLong=False) # only ingroup; extended func now also works with visits;
    outchimerasplancnti_gp = outchimerasplancnti[goodpos2useTree] # ATCGN translation; outchimerasplancnti is chimera of patient-specific outgroup and goodpos that remained unresolved (ie. 4) are substituted by major allele of all outgroup samples.
    numBaseGenome_covThresholdIngroup = numBaseGenome_covThreshold[ingroup_bool]
    patient_time_since_v1 = pat_visit_timing[pat_visit_timing['patient']==int(subject_fld_label.split('_')[2])][['time1','time2','time3','time4','time5']].values.tolist()[0] # adds actual timing of visit (not always right on time 0/1/2/3/9 months)

    # loop over each sample and assess #diffs to outgroup
    diff_count = np.zeros(visits[ingroup_bool].shape,dtype=float64)
    for array in list_visit_idx:
        for sample_idx in array:
            counter = 0
            for i,out_nt in enumerate( outchimerasplancnti_gp ):
                if out_nt != 4 and calls_for_treei_ingroup[i,sample_idx] != 4 and out_nt != calls_for_treei_ingroup[i,sample_idx]:
                    counter += 1 # counts nt distance sample and ancnt
            diff_count[sample_idx] = counter

    # get counts 
    hdr = np.array(['Month','rate','numBaseGenome'])
    data = np.zeros( shape=(len(diff_count),3),dtype=float)
    mean_count_per_visit = []
    for i,v in enumerate(list_visit_idx):
        if v.size != 0:
            data[v,1] = diff_count[v]
            data[v,0] = patient_time_since_v1[i] #[0,1,2,3,9][i] # time-in-months corresponds to #visit
            data[v,2] = numBaseGenome_covThresholdIngroup.iloc[v].values[:,0]
            mean_count_per_visit.append( np.array([ patient_time_since_v1[i],mean(data[v,1]/data[v,2]) ]) )
    mean_count_per_visit = np.array(mean_count_per_visit)
    num_samples = data.shape[0]
    num_visits = len(np.unique(data[:,0]))
    
    # save data
    mol_clock_data_dc[subject_fld_label] = data # store data for subsequent regression analysis that includes correction for genome-wide base calls


    # %% PLOT seaborn
    # calc regression
    # Warning: make conditional for data with > 1 visit! adjust plotting function, too.
    regress_dict = apy.calc_regression(np.column_stack((data[:,0],data[:,1]/data[:,2])),subject_fld_label,num_samples,num_visits)
    regression_res_lod.append( regress_dict ) # store, turn to pd dataframe for all subjects
    
    basescale = 1000000 # rescale genome-wide molecular rate to human-readable numbers
    apy.plot_molecular_clock(np.column_stack( (data[:,0],(data[:,1]/data[:,2]) )),hdr[0:2],mean_count_per_visit,regress_dict,basescale,subject_fld_label) # plot in pdf/molclock/
    
    # %% 
    # =============================================================================
    #     SNP frequency change Module
    # =============================================================================
    # record all SNPs w/ AF change > X; see "allele_freq_change_min"
    # NOTE: AF change based on difference to outgroup_spl allele (recorded anc in annotation_mutation). if out_sample == 4, use ref_nt. 
    allele_freq_change_min = 0.3 # the minimum shift in observed allele freq between visits required for reporting
    for idx,anno_gene in annotation_mutations.iterrows():
        freq_visit = np.zeros(5)
        count_spl_visit = np.zeros(5,dtype=int)
        ref_nti = np.where(NTs == anno_gene['nt_anc'])[0]
        if ref_nti.size == 0: # if 'nt_anc' == 4, use ref_nt. alternative is to skip site (but I rather re-polarize it if needed than missing a potential interesting mutation).
            ref_nti = np.where(NTs == anno_gene['nt_ref'])[0]
        alleles_i = calls[ goodpos2use[idx] , ingroup_bool ]
        subj_visits = visits_all[ goodsamples ] # aka visits prev defined
        subj_visits = subj_visits[ ingroup_bool ] # remove outgroup samples
        for v in sort(np.unique(subj_visits)): # loop through all visits available
            allelesi_at_visits = alleles_i[ subj_visits==v ] 
            if sum(allelesi_at_visits != 4) != 0: # test denominator not 0 (all samples 4)
                freq_visit[v-1] = np.sum(np.all(((allelesi_at_visits != 4) ,(allelesi_at_visits != ref_nti)),axis=0)) / np.sum(allelesi_at_visits != 4) # visit 1-based. Array 0-based
            count_spl_visit[v-1] = sum(subj_visits==v)
        # calc max allele frequency chronological (min >> max from vi to vn; subjects with only one visit not considered)
        bool_valid_visit = count_spl_visit > 2
        if sum(bool_valid_visit) > 1:
            valid_freq_visit = freq_visit[ bool_valid_visit ]
            ## initially removed mut that had peak at v1
            max_i = np.argmax(valid_freq_visit)
            # if max_i != 0: # if max already at first valid visit >> skip                
            #     freq_min = min(valid_freq_visit[:max_i]) # get min from values prior max freq visit!                
            freq_min = min(valid_freq_visit) # get min from values prior max freq visit!                    
            freq_max = valid_freq_visit[max_i]
            af_change = freq_max - freq_min
            if af_change > allele_freq_change_min:
                cand_mut_dict = {}
                # build dict for pd.dataframe
                cand_mut_dict['0chr'] = anno_gene['chr']
                cand_mut_dict['0pos'] = anno_gene['pos']
                cand_mut_dict['0tag'] = anno_gene['locustag']
                cand_mut_dict['0visit1_freqs'] = freq_visit[0]
                cand_mut_dict['0visit2_freqs'] = freq_visit[1]
                cand_mut_dict['0visit3_freqs'] = freq_visit[2]
                cand_mut_dict['0visit4_freqs'] = freq_visit[3]
                cand_mut_dict['0visit5_freqs'] = freq_visit[4]
                cand_mut_dict['0visit1_counts'] = count_spl_visit[0]
                cand_mut_dict['0visit2_counts'] = count_spl_visit[1]
                cand_mut_dict['0visit3_counts'] = count_spl_visit[2]
                cand_mut_dict['0visit4_counts'] = count_spl_visit[3]
                cand_mut_dict['0visit5_counts'] = count_spl_visit[4]
                cand_mut_dict['0subject'] = subject_fld_label
                for col in annotation_mutations.columns: # get all data from dataframe
                    cand_mut_dict[col] = anno_gene[col]
                snp_freq_shifts_anno_lod.append(cand_mut_dict)
    
    ## change back to base folder
    os.chdir('../..')

    # %%
######
######
###### END OF LOOP!
######
######
          

# assumes located in base folder right now 'within_patient_analysis'
mean_cov_p_per_sample = np.array(['sample','mean_p_cov'])


# %% Write table with annotation_mutation for all para evo candidates
all_para_evo_mut = pd.DataFrame() # Create an empty dataframe, this will be your final dataframe
for key, sub_df in annotation_mutation_allParaSignal.items():
    new_df_wPatID = sub_df.assign(patient=key)
    all_para_evo_mut = all_para_evo_mut.append(new_df_wPatID, ignore_index=False,sort=True) # Add your sub_df one by one
all_para_evo_mut.to_csv('parallelEvo_cand_annotation_mutation.csv')

# %% write mol clock data as json
with open('molecular_clock_data.pk1', 'wb') as fp:
    pickle.dump(mol_clock_data_dc, fp)


# =============================================================================
# %% Infer presence and mutation type of cand genes across all subjects
# =============================================================================
para_cand_allS_type = pd.DataFrame(columns=mut_genes_df.keys(),index=np.unique(np.array([x for x in para_evo_cand if not pd.isnull(x)])) ) # geeky np.unique required bcs numpy cannot handle mix string and nan

# loop over each gene
for cand in para_cand_allS_type.index.values:
    for pat in para_cand_allS_type.columns.values:
        subset = mut_genes_df[pat].loc[mut_genes_df[pat]['orthologtag'] == cand]
        if len(subset.index.values) > 0:
            para_cand_allS_type.at[cand,pat] = "".join(subset['type'].values)
        elif 'orthologtag1' in mut_genes_df[pat].columns.values: # orthologtag1/2 only reported in set at least for one SNP
            subset = mut_genes_df[pat].loc[mut_genes_df[pat]['locustag1'] == cand]
            if len(subset.index.values) > 0:
                para_cand_allS_type.at[cand,pat] = "".join(subset['type'].values)
            else:
                subset = mut_genes_df[pat].loc[mut_genes_df[pat]['orthologtag2'] == cand]
                if len(subset.index.values) > 0:
                    para_cand_allS_type.at[cand,pat] = "".join(subset['type'].values)
para_cand_allS_type = para_cand_allS_type.fillna('-')  # replace NaN with "-" to increase visibility
para_cand_allS_type.to_csv('parallelEvo_cand_allSubjects.csv')

# get gene annotation info for candidates from mut_genes_df
mut_cand_subjects_ls = []
for tag in para_cand_allS_type.index.values:
    for k in mut_genes_df.keys():
        if (mut_genes_df[k]['orthologtag'] == tag ).any():
            row_df = mut_genes_df[k][ mut_genes_df[k]['orthologtag'] == tag ]
            row_df.insert(2,"ssdubject",[k]*row_df.shape[0]  ) # add subject column; *row_df.shape[0] repeat value up to nrow(df_row)
            row_df.insert(2,"genelength", list(row_df.loc[:,'loc2'] - row_df.loc[:,'loc1']) ) # add subject column; 
            mut_cand_subjects_ls.append(row_df)
mut_cand_subjects_df = pd.concat(mut_cand_subjects_ls)

# ortholog info for candidates
ortholog_df_cand = ortholog_df.loc[ para_cand_allS_type.index.values , : ]
ortholog_df_cand.to_csv('parallelEvo_ortholog_cand_allSubjects.csv')

# save all mutations in csv (ie. all goodpos); add subject identifier first to df
for subj in mut_genes_df.keys():
    mut_genes_df[subj]['subject'] = subj


pd.concat(mut_genes_df,axis=0, ignore_index=True, sort=True).to_csv('allGoodpos_allSubjects.csv') # sort=True to align colmns which are not always the same, dependeing if genic/nongenic SNPs part of goodPos per subject

# %% Write table with molecular clock data
df_reg = pd.DataFrame(regression_res_lod)
df_reg.to_csv('mol_clock_regression_results.csv')



# 
# =============================================================================
# %% Analyse SNP frequency record data
# =============================================================================
snp_freq_shifts_anno = pd.DataFrame(snp_freq_shifts_anno_lod) 

# for I/P mutations: turn all to I (P somewhat random); add up/downstream info to 0tag and gene description to product
visit_counts = ['0visit1_counts','0visit2_counts','0visit3_counts','0visit4_counts','0visit5_counts']
visit_freqs = ['0visit1_freqs','0visit2_freqs','0visit3_freqs','0visit4_freqs','0visit5_freqs']
rising_mut_only = []

## generate idx file containing only mut that rise in freq during sampling
for i,row in snp_freq_shifts_anno.iterrows():
    if row.isnull()['0tag']: # snps outside of annotated gene
        if not row.isnull()['locustag1'] and not row.isnull()['locustag2']: # only process snps that have gene up or downstream [removes hits on small contigs w/o annotation]
            snp_freq_shifts_anno.at[i,'0tag'] = str(int(row['distance1']))+":"+row['locustag1']+";"+str(int(row['distance2']))+":"+row['locustag2']
            snp_freq_shifts_anno.at[i,'gene'] = row['gene1']+";"+row['gene2']
            snp_freq_shifts_anno.at[i,'product'] = row['product1']+";"+row['product2']
            snp_freq_shifts_anno.at[i,'type'] = 'I'
        else:
            print(i) # SNP located on contig w/o any annotation
    # filter SNPs for max freq at last visit w/ data (ie. only contains mut that rise in freq!)
    visits_w_data_bool = (row[visit_counts] > 0).values
    visit_w_data_freq = row[visit_freqs][visits_w_data_bool].values
    if np.max(visit_w_data_freq)==visit_w_data_freq[-1]:
        rising_mut_only.append(i)
                     
snp_freq_shifts_anno.to_csv('freqEvo_cand_allSubjects.csv')
snp_freq_shifts_anno=snp_freq_shifts_anno.iloc[rising_mut_only]
snp_freq_shifts_anno.to_csv('freqEvo_cand_allSubjects_risingFreqMutOnly.csv')

# Obtain genes with mutational sweep in multiple patients
orthologs,ortholog_sweep_count = np.unique(snp_freq_shifts_anno['orthologtag'].dropna().values, return_counts=True)
orthologs[ortholog_sweep_count>1]

## check for cand from para_cand_allS_type
cand_freq_and_para_signature = snp_freq_shifts_anno.loc[ snp_freq_shifts_anno['orthologtag'].isin(para_cand_allS_type.index),: ]
cand_freq_and_para_signature.to_csv('freqEvo_paraEvo_shared_cand.csv')


# =============================================================================
# %% Write all para evo candidates 
# =============================================================================

with open('all_patients_parallel_evo_results.tsv','w+') as file:
    for s in allP_gene_paraevo.keys():
        allP_gene_paraevo[s].to_csv(file,sep="\t",index=False,header=True)

