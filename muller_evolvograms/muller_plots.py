#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 12:21:42 2021

@author: u_key
"""

############################################################################################################################################
############################################################################################################################################
########## Evolvograms / muller plots
############################################################################################################################################
############################################################################################################################################
# =============================================================================
# %% Available files
# =============================================================================
# muller_conda_environment.yml
# ../metadata/patient_visit_date_timepast.csv
# freqEvo_cand_allSubjects.csv # provided, but otherwise generated with within_patient_analysis.py

## In addition, within each subject folder different input and intermediate files are provided.

# =============================================================================
# %% Setup
# =============================================================================
# use the conda env provided. alternatively, generate environemnt incl. all dependencies urself using the recipe below:

# use the pypi muller package
# https://github.com/cdeitrick/Lolipop


## The installation of pygraphviz (and its dependecies) is somewhat in conflict with conda. Thus I circumvent it via brew.
## publish conda env and add single non-conda steps
# ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
# conda create -n muller python=3.7
# conda activate muller
# brew install graphviz  
# conda install scipy loguru matplotlib pandas seaborn  xlrd --channel conda-forge --channel bioconda
# pip install lolipop # v0.9.0
# conda install -c r r r-ggplot2 
# conda update -c rdonnellyr -c main --all # R package failed to load. Found this fix: https://github.com/conda/conda/issues/6183

## in R
# install.packages("ggmuller") # https://github.com/robjohnnoble/ggmuller


# =============================================================================
# %% Modules
# =============================================================================

import sys,os,re,gzip,json
import glob,pickle,subprocess,random,time,datetime,collections,warnings
import numpy as np
from scipy import stats
import pandas as pd

# =============================================================================
# %% Functions
# =============================================================================


def smooth_nested_muller_freq_table(arr,hdr,edges,traj_genotype):
    ''' Simulate exponential frequency shifts for muller plots '''
    # returns new (possibly extended) array and array with new generation ids (colnames for muller.table)
    # arr & traj_genotype need to have same vertical order!
    # each time point (visit) subdivided into 100 steps. 5 nested level each pushed by 20 steps. 40 steps to go from start to end freq between two visits
  
    nested_level_rise_fix = [[(0,100)],[(0,80),(20,100)],[(0,60),(20,80),(40,100)],[(0,40),(20,60),(40,80),(60,100)],[(0,20),(20,40),(40,60),(60,80),(80,100)]] # lol of tuples with start end rise o fixation
    nested_level_rise_itm = [[0],(0,50),(0,33,66),(0,25,50,75),(0,20,40,60,80)] # list of tuples that define push of freq rise for mut that go until intermediate freq; element 0 is []bcs py does drop tuple when only single item
    
    ncol_arr = arr.shape[1]
    g = 100 # number of generations we extend each time step
    ## get nested level for each mut in arr
    nested_level_traj = np.ones(traj_genotype.shape[0]) # record nested level of genotype from 1 to x. same order as arr/traj_genotype
    for j,row in enumerate(traj_genotype):
        # loop over each traj/mut and infer nest level 
        nest = traj_genotype[j,1] # get genotype id
        parent = edges[ edges[:,1]==nest , 0 ][0] # get parental genotype id
        while parent != 'genotype-0':
            # traverse through edges until genotype-0 (origin) to infer nest level
            nested_level_traj[j] = nested_level_traj[j] + 1
            nest = parent # get genotype id
            parent = edges[ edges[:,1]==nest , 0 ] # get parental genotype id
    nested_level_traj = nested_level_traj.astype(int)
    ## add simulations
    for i in np.arange((ncol_arr-1)): # loop over each but last col/time point        
        if not np.equal(arr[:,i],arr[:,(i+1)]).all():
            # if any value between two time points does not agree
            ## infer which mut change type, necessary for correct sim
            rise_fix_bool = np.logical_and(np.invert(np.equal(arr[:,i],arr[:,(i+1)])),arr[:,(i+1)] == 1) # rise to fixation
            rise_im_bool = np.invert(np.equal(arr[:,i],arr[:,(i+1)]))*(arr[:,(i)] < arr[:,(i+1)])*(arr[:,(i+1)]!=1) # a*b*c is 'and' for boolean of same length!, rise to im freq
            decline_bool = np.logical_and(np.invert(np.equal(arr[:,i],arr[:,(i+1)])),arr[:,(i)] > arr[:,(i+1)]) # decline
            mut_type = np.repeat('flat',arr.shape[0])
            mut_type[rise_fix_bool] = 'fix' # rise_fix
            mut_type[rise_im_bool] = 'itm' # rise to intermediate
            mut_type[decline_bool] = 'decl' # decline
            extra_arr = np.array([]) # array to collect sims
            for j in np.arange(arr.shape[0]):
                # loop over each traj/mut and simulate exponential curve accordingly
                print(i,j)
                # get freqs
                f_s = arr[j,i] # frequency start
                f_e = arr[j,(i+1)] # frequency end
                # freq shift of preexisting variation: nested level defines freq reach
                if f_e == 0:
                    f_e = 0.001 # ln error for 0
                if f_s == 0:
                    f_s = 0.001 # ln error for 0
                if mut_type[j] == 'flat':
                    # no change. fill gen
                    ynew = np.repeat(f_s,g)
                elif mut_type[j] == 'fix':
                    ## mut goes to fixation: nestedness affect start/end of freq shift. start @ 0 if start_freq != 0
                    mut_nest_lvl = nested_level_traj[j]
                    min_nest_level_iteration = np.min(np.unique(nested_level_traj[rise_fix_bool])) # get the lowest nest level that changes between current time point i and i+1
                    mut_nest_lvl = mut_nest_lvl - min_nest_level_iteration # correct nest level by relevant nest levels in iteration. also 0-base mut_nest_level!
                    num_nest_levels_iteration = np.unique(nested_level_traj[rise_fix_bool]).size - 1 # number of nest levels for rise_fix. 0-based
                    start_g,end_g = nested_level_rise_fix[num_nest_levels_iteration][mut_nest_lvl]
                    if f_s > 0.001:
                        start_g = 0 # start freq rise at 0
                    num_g_sim = end_g - start_g
                    x = np.exp(((np.log(f_e)-np.log(f_s))/num_g_sim)) # x = e^((ln(f_e) - ln(f_s) ) / g)
                    ynew = np.concatenate(( np.repeat(f_s,start_g) , np.array([f_s*x**y for y in np.arange(1,num_g_sim+1) ]), np.repeat(f_e,g-end_g) ))
                elif mut_type[j] == 'itm':
                    ## mut rises to intermediate freq. start of rise is nested, but not end
                    mut_nest_lvl = nested_level_traj[j]
                    min_nest_level_iteration = np.min(np.unique(nested_level_traj[rise_im_bool])) # get the lowest nest level that changes between current time point i and i+1 for type_mut
                    mut_nest_lvl = mut_nest_lvl - min_nest_level_iteration # correct nest level by relevant nest levels in iteration. also 0-base mut_nest_level!
                    num_nest_levels_iteration = np.unique(nested_level_traj[rise_fix_bool]).size - 1 # number of nest levels for rise_fix. 0-based
                    if np.any(arr[rise_fix_bool,i]==0):
                        ## correct for fixed mutations. Nesting needs to consider it (S9) but only f fixed mut from 0 to 1!
                        num_nest_level_riseFix = np.unique(nested_level_traj[rise_fix_bool]).size - 1 # num genotypes that get fixed. 0-based
                    else:
                        num_nest_level_riseFix = 0
                    start_g = nested_level_rise_itm[num_nest_levels_iteration+num_nest_level_riseFix][mut_nest_lvl+num_nest_level_riseFix]
                    num_g_sim = g - start_g
                    x = np.exp(((np.log(f_e)-np.log(f_s))/num_g_sim)) # x = e^((ln(f_e) - ln(f_s) ) / g)
                    ynew = np.concatenate(( np.repeat(f_s,start_g) , np.array([f_s*x**y for y in np.arange(1,num_g_sim+1) ]) ))
                elif mut_type[j] == 'decl':
                    ## mut declines. decline finished when first rise2fixation finished.
                    # get number of nested gt that get fixed; decline horizontal push based on #nested levels getting fixed
                    if np.any(rise_fix_bool):
                        ## there is at least one gt going to fixation. decline will be lienar
                        num_nest_level_riseFix = np.unique(nested_level_traj[rise_fix_bool]).size - 1 # num genotypes that get fixed. 0-based
                        end_g = nested_level_rise_fix[num_nest_level_riseFix][0][1] # endpoint of mutation fixed first in that iteration
                        num_g_sim = end_g
                        # x = np.exp(((np.log(f_e)-np.log(f_s))/num_g_sim)) # x = e^((ln(f_e) - ln(f_s) ) / g)
                        # ynew = np.concatenate(( np.flip(-1*np.array([f_s*x**y for y in np.arange(1,num_g_sim+1) ]))+f_s+f_e , np.repeat(f_e,g-end_g) ))
                        ynew = np.concatenate(( np.linspace(f_s,f_e,num_g_sim) , np.repeat(f_e,g-end_g) ))
                    else:
                        ## there is no gt going to fixation. decline will be exponential
                        num_g_sim = g
                        # x = np.exp(((np.log(f_e)-np.log(f_s))/num_g_sim)) # x = e^((ln(f_e) - ln(f_s) ) / g)
                        x = np.exp(((np.log(f_e)-np.log(f_s))/num_g_sim)) # x = e^((ln(f_e) - ln(f_s) ) / g)
                        ynew = np.flip(-1*np.array([f_s*x**y for y in np.arange(1,num_g_sim+1) ]))+f_s+f_e
                        # ynew = np.array([f_s*x**y for y in np.arange(1,num_g_sim+1) ])
                ## collect sim data   
                extra_arr = np.concatenate((extra_arr,ynew),axis=0) # growth/decline simulated data.   
            ## build extended array and hdr that incorporate simulations 
            extra_arr = extra_arr.reshape(arr.shape[0],ynew.size)
            extra_hdr = np.round(np.linspace(start=hdr[i],stop=hdr[i+1],num=ynew.size+1 ),0)[:-1] # adjust float to integer values (require integer generation time); use linspace to be flexible with 100 or 600gen (v4>v5)
        else:
            extra_arr = arr[:,i:(i+1)]
            extra_hdr = np.array([hdr[i]])
        if i == 0:
            res = extra_arr
            new_hdr = extra_hdr
        else:
            res = np.hstack([res,extra_arr])
            new_hdr = np.append(new_hdr,extra_hdr)
    res = np.hstack([res,arr[ :,(ncol_arr-1):ncol_arr ]]) # append last column arr
    new_hdr = np.append(new_hdr,hdr[-1]) # append last col name (gen time)
    return pd.DataFrame(data=res,columns=(new_hdr.astype(int)))

# =============================================================================
# %% Working directory and data
# =============================================================================
## read SNP freq table
data = pd.read_csv('freqEvo_cand_allSubjects.csv') # all SNVs with freq diff >0.3

## read visit times
visit_time = pd.read_csv('../metadata/patient_visit_date_timepast.csv')

# %%
##############################
###### subject_staphAD_4
##############################

subj_id = "subject_staphAD_4"
visit_time_pat = visit_time.loc[visit_time['patient'] == int(subj_id.split('_')[2])][['time1','time2','time3','time4','time5']].round(0).astype(int)


## output folder
res_folder = subj_id
os.chdir(res_folder)

## get muller/lolipop input table
lod_freq = []
ctr = 1
for index, row in data.iterrows():
    if row['0subject'] == subj_id:
        dd = {}
        dd['Population'] = row['0subject']
        dd['Trajectory'] = ctr
        dd['Chromosome'] = row['chr']
        dd['Position'] = row['pos']
        print(row['pos'])
        dd['Class'] = 'SNP'
        dd['Mutation'] = str(row['chr'])+"_"+str(row['pos'])+"_"+str(row['muts'])
        for j,i in enumerate(range(1,6)): # loop through visit freq data
            # print(i)
            if row['0visit'+str(i)+'_counts'] != 0: # only record visit when samples avail (V1 always has samples!); 
                dd[str(visit_time_pat.values[0,j])] = row['0visit'+str(i)+'_freqs']
        if dd['Position'] != 43772: # 43772 position that has no SNP in visits with > 3 samples. terrible fix
            ctr = ctr + 1
            lod_freq.append(dd)
muller_tbl = pd.DataFrame(lod_freq)

## filter SNP freqs
muller_tbl = muller_tbl.drop(['3'], axis=1) # REMOVE visit 4, with only 1 isolate
## Due to a bug in lolipop V5 emerging clades need to be adjusted for correct repesentation in muller; in order to represent both emerging strains at v5 they have to be set to 0.9 and 1, respectively (despite extensive debugging I found no alternative)
# note: the trajectories need to be adjusted for individual frequency tables!
muller_tbl['12'][[0,6,7,8]]=0.9
muller_tbl['12'][[3,4,5]]=1
muller_tbl.to_csv("muller_tbl.tsv", sep="\t",index=False)



## run lolipop
# lolipop does not always find correct ancestry and genotypes. Perform iterative process to fix by manually building genotypes and ancestry input files
subprocess.run('rm -r results'.split(' '))
subprocess.run('lolipop lineage --input muller_tbl.tsv -o results'.split(' '))
subprocess.run('lolipop lineage --input muller_tbl.tsv -o results --known-ancestry known_ancestry.tsv --known-genotypes known_genotypes.csv'.split(' '))


## smooth out frequency rise by exponential function
# turn 0,1,2,3,9 > 0,100,200,300,900
# every rise in freq preceded by 10 step exponential cline (drop)
freq_arr = np.array(muller_tbl.iloc[:,6:12])
arr = np.round(freq_arr,2)
hdr = muller_tbl.columns[6:12].astype('int') * 100


## read in trajectory > genotype array and sort by trajectory numeric id...ie. sort as mueller_tbl [ony available after first run of lolipop]
with open('results/supplementary-files/muller_tbl.genotypemembers.json') as f:
  genotype_traj_dc = json.load(f)

traj_genotype = np.array([])
for gt in genotype_traj_dc:
    for t in genotype_traj_dc[gt]:
        traj_genotype = np.append(traj_genotype,[int(t),gt])
traj_genotype = traj_genotype.reshape(int(traj_genotype.size/2),2)
traj_genotype=traj_genotype[ np.argsort(traj_genotype[:,0].astype('int')) ,:] # sort based on traj value, which turns array in same order as freq_tbl/muller_tbl


## read in edges [ony available after first run of lolipop]
edges = np.loadtxt('results/tables/muller_tbl.edges.tsv',skiprows=1,delimiter='\t',dtype=object)

## smooth freq shifts and consider nest-level for lolipop                
extended_data_df = smooth_nested_muller_freq_table(arr,hdr,edges,traj_genotype)
# manually changed genotype 4 @ V5: 0.9 and genotype-2 @ V5: 1

extended_muller_tbl = pd.concat([muller_tbl.iloc[:,0:6], extended_data_df], axis=1, sort=False)
extended_muller_tbl.to_csv("muller_extended_tbl.tsv", sep="\t",index=False)

##run lolipop from cl
subprocess.run('rm -r results_nested'.split(' '))
subprocess.run('lolipop lineage --input muller_extended_tbl.tsv -o results_nested --detection 0.00 --significant 0.00 --known-genotypes known_genotypes.csv --known-ancestry known_ancestry.tsv'.split(' '))

## polish muller plot using ggplot/ggmuller in R
subprocess.run('Rscript plot_polished_muller_p4.r')
# > generates figure (pdf) in working directory
                
# leave working subject-specific directory
os.chdir('..')

##############################
###### subject_staphAD_9
##############################

subj_id = "subject_staphAD_9"
visit_time_pat = visit_time.loc[visit_time['patient'] == int(subj_id.split('_')[2])][['time1','time2','time3','time4','time5']].round(0).astype(int)


## output folder
res_folder = subj_id
os.chdir(res_folder)


## get muller/lolipop input table
lod_freq = []
ctr = 1
for index, row in data.iterrows():
    if row['0subject'] == subj_id:
        dd = {}
        dd['Population'] = row['0subject']
        dd['Trajectory'] = ctr
        dd['Chromosome'] = row['chr']
        dd['Position'] = row['pos']
        print(row['pos'])
        dd['Class'] = 'SNP'
        dd['Mutation'] = str(row['chr'])+"_"+str(row['pos'])+"_"+str(row['muts'])
        for j,i in enumerate(range(1,6)): # loop through visit freq data
            if row['0visit'+str(i)+'_counts'] != 0: # only record visit when samples avail (V1 always has samples!); 
                dd[str(visit_time_pat.values[0,j])] = row['0visit'+str(i)+'_freqs']
        ctr = ctr + 1
        lod_freq.append(dd)
muller_tbl = pd.DataFrame(lod_freq)
muller_tbl.reset_index(drop=True, inplace=True)
muller_tbl.to_csv("muller_tbl.tsv", sep="\t",index=False)


## run lolipop in conda base
# fix genotypes and ancestry if needed and append to input files 
subprocess.run('rm -r results'.split(' '))
subprocess.run('lolipop lineage --input muller_tbl.tsv -o results --known-genotypes known_genotypes.csv --known-ancestry known_ancestry.tsv'.split(' '))


## smooth out frequency rise by exponential function
# turn 0,1,2,3,9 > 0,100,200,300,900
# every rise in freq preceded by 10 step exponential cline (drop)

freq_arr = np.array(muller_tbl.iloc[:,6:12])
arr = freq_arr
hdr = muller_tbl.columns[6:12].astype('int') * 100


## read in trajectory > genotype array and sort by trajectory numeric id...ie. sort as mueller_tbl [ony available after first run of lolipop]
with open("muller/"+subj_id+'/results/supplementary-files/muller_tbl.genotypemembers.json') as f:
  genotype_traj_dc = json.load(f)

traj_genotype = np.array([])
for gt in genotype_traj_dc:
    for t in genotype_traj_dc[gt]:
        traj_genotype = np.append(traj_genotype,[int(t),gt])
traj_genotype = traj_genotype.reshape(int(traj_genotype.size/2),2)
traj_genotype=traj_genotype[ np.argsort(traj_genotype[:,0].astype('int')) ,:] # sort based on traj value, which turns array in same order as freq_tbl/muller_tbl


## read in edges [ony available after first run of lolipop]
edges = np.loadtxt('muller/'+subj_id+'/results/tables/muller_tbl.edges.tsv',skiprows=1,delimiter='\t',dtype=object)


## smooth freq shifts and consider nest-level for lolipop                
extended_data_df = smooth_nested_muller_freq_table(arr,hdr,edges,traj_genotype)


extended_muller_tbl = pd.concat([muller_tbl.iloc[:,0:6], extended_data_df], axis=1, sort=False,ignore_index=True)
extended_muller_tbl.to_csv("muller/"+subj_id+"/muller_extended_tbl.tsv", sep="\t",index=False)

##run lolipop from cl
subprocess.run('rm -r results_nested'.split(' '))
subprocess.run('lolipop lineage --input muller_extended_tbl.tsv -o results_nested --detection 0.00 --significant 0.00 --known-genotypes known_genotypes.csv --known-ancestry known_ancestry.tsv'.split(' '))


## polish muller plot using ggplot/ggmuller in R
subprocess.run('Rscript plot_polished_muller_p4.r')
# > generates figure (pdf) in working directory

# leave working subject-specific directory
os.chdir('..')

##############################
###### subject_staphAD_12
##############################

subj_id = "subject_staphAD_12"
visit_time_pat = visit_time.loc[visit_time['patient'] == int(subj_id.split('_')[2])][['time1','time2','time3','time4','time5']].round(0).astype(int)

## output folder
res_folder = subj_id
os.chdir(res_folder)

## get muller/lolipop input table
lod_freq = []
ctr = 1
for index, row in data.iterrows():
    if row['0subject'] == subj_id:
        dd = {}
        dd['Population'] = row['0subject']
        dd['Trajectory'] = ctr
        dd['Chromosome'] = row['chr']
        dd['Position'] = row['pos']
        print(row['pos'])
        dd['Class'] = 'SNP'
        dd['Mutation'] = str(row['chr'])+"_"+str(row['pos'])+"_"+str(row['muts'])
        for j,i in enumerate(range(1,6)): # loop through visit freq data
           if row['0visit'+str(i)+'_counts'] != 0: # only record visit when samples avail (V1 always has samples!); 
                dd[str(visit_time_pat.values[0,j])] = row['0visit'+str(i)+'_freqs']
        ctr = ctr + 1
        lod_freq.append(dd)
muller_tbl = pd.DataFrame(lod_freq)
muller_tbl.to_csv("muller_tbl.tsv", sep="\t",index=False)



## run lolipop in conda base
# fix genotypes and ancestry if needed and append to input files
subprocess.run('rm -r results'.split(' '))
subprocess.run('lolipop lineage --input muller_tbl.tsv -o results --known-genotypes known_genotypes.csv --known-ancestry known_ancestry.tsv'.split(' '))


## smooth out frequency rise by exponential function
# turn 0,1,2,3,9 > 0,100,200,300,900
# every rise in freq preceded by 10 step exponential cline (drop)

freq_arr = np.array(muller_tbl.iloc[:,6:12])
arr = np.round(freq_arr,2)
hdr = muller_tbl.columns[6:12].astype('int') * 100


## read in trajectory > genotype array and sort by trajectory numeric id...ie. sort as mueller_tbl [ony available after first run of lolipop]
with open('results/supplementary-files/muller_tbl.genotypemembers.json') as f:
  genotype_traj_dc = json.load(f)

traj_genotype = np.array([])
for gt in genotype_traj_dc:
    for t in genotype_traj_dc[gt]:
        traj_genotype = np.append(traj_genotype,[int(t),gt])
traj_genotype = traj_genotype.reshape(int(traj_genotype.size/2),2)
traj_genotype=traj_genotype[ np.argsort(traj_genotype[:,0].astype('int')) ,:] # sort based on traj value, which turns array in same order as freq_tbl/muller_tbl


## read in edges [ony available after first run of lolipop]
edges = np.loadtxt('results/tables/muller_tbl.edges.tsv',skiprows=1,delimiter='\t',dtype=object)

## smooth freq shifts and consider nest-level for lolipop                
extended_data_df = smooth_nested_muller_freq_table(arr,hdr,edges,traj_genotype)
extended_data_df[extended_data_df > 1] = 1

# extended_data_norm_df = extended_data_df.div(extended_data_df.sum(axis=0), axis=1)

extended_muller_tbl = pd.concat([muller_tbl.iloc[:,0:6], extended_data_df], axis=1, sort=False)
extended_muller_tbl.to_csv("muller_extended_tbl.tsv", sep="\t",index=False)

##run lolipop from cl
subprocess.run('rm -r results_nested'.split(' '))
subprocess.run('lolipop lineage --input muller_extended_tbl.tsv -o results_nested --known-ancestry known_ancestry.tsv --detection 0.00 --significant 0.00 --known-genotypes known_genotypes.csv'.split(' '))

## polish muller plot using ggplot/ggmuller in R
subprocess.run('Rscript plot_polished_muller_p12.r')
# > generates figure (pdf) in working directory

# leave working subject-specific directory
os.chdir('..')

##############################
###### subject_staphAD_15
##############################

subj_id = "StaphAureus_fmk_15"
visit_time_pat = visit_time.loc[visit_time['patient'] == int(subj_id.split('_')[2])][['time1','time2','time3','time4','time5']].round(0).astype(int)

## output folder
res_folder = subj_id
os.chdir(res_folder)

## get muller/lolipop input table
lod_freq = []
ctr = 1
for index, row in data.iterrows():
    if row['0subject'] == subj_id:
        dd = {}
        dd['Population'] = row['0subject']
        dd['Trajectory'] = ctr
        dd['Chromosome'] = row['chr']
        dd['Position'] = row['pos']
        print(row['pos'])
        dd['Class'] = 'SNP'
        dd['Mutation'] = str(row['chr'])+"_"+str(row['pos'])+"_"+str(row['muts'])
        for j,i in enumerate(range(1,6)): # loop through visit freq data
           if row['0visit'+str(i)+'_counts'] != 0: # only record visit when samples avail (V1 always has samples!); 
                dd[str(visit_time_pat.values[0,j])] = row['0visit'+str(i)+'_freqs']
        ctr = ctr + 1
        lod_freq.append(dd)
muller_tbl = pd.DataFrame(lod_freq)
muller_tbl = muller_tbl.iloc[[1,4]]
muller_tbl.reset_index(drop=True, inplace=True)
muller_tbl = muller_tbl.drop(['1'], axis=1) # REMOVE visit 2, with only 2 isolate
muller_tbl.to_csv("muller_tbl.tsv", sep="\t",index=False)



## run lolipop in conda base
# fix genotypes and ancestry if needed and append to input files
subprocess.run('rm -r results'.split(' '))
subprocess.run('lolipop lineage --input muller_tbl.tsv -o results'.split(' '))


## smooth out frequency rise by exponential function
# turn 0,1,2,3,9 > 0,100,200,300,900
# every rise in freq preceded by 10 step exponential cline (drop)

freq_arr = np.array(muller_tbl.iloc[:,6:12])
arr = np.round(freq_arr,2)
hdr = muller_tbl.columns[6:12].astype('int') * 100


## read in trajectory > genotype array and sort by trajectory numeric id...ie. sort as mueller_tbl [ony available after first run of lolipop]
with open('results/supplementary-files/muller_tbl.genotypemembers.json') as f:
  genotype_traj_dc = json.load(f)

traj_genotype = np.array([])
for gt in genotype_traj_dc:
    for t in genotype_traj_dc[gt]:
        traj_genotype = np.append(traj_genotype,[int(t),gt])
traj_genotype = traj_genotype.reshape(int(traj_genotype.size/2),2)
traj_genotype=traj_genotype[ np.argsort(traj_genotype[:,0].astype('int')) ,:] # sort based on traj value, which turns array in same order as freq_tbl/muller_tbl


## read in edges [ony available after first run of lolipop]
edges = np.loadtxt('results/tables/muller_tbl.edges.tsv',skiprows=1,delimiter='\t',dtype=object)

## smooth freq shifts and consider nest-level for lolipop                
extended_data_df = smooth_nested_muller_freq_table(arr,hdr,edges,traj_genotype)
extended_data_df[extended_data_df > 1] = 1

# extended_data_norm_df = extended_data_df.div(extended_data_df.sum(axis=0), axis=1)

extended_muller_tbl = pd.concat([muller_tbl.iloc[:,0:6], extended_data_df], axis=1, sort=False)
extended_muller_tbl.to_csv("muller_extended_tbl.tsv", sep="\t",index=False)

##run lolipop from cl
subprocess.run('rm -r results_nested'.split(' '))
subprocess.run('lolipop lineage --input muller_extended_tbl.tsv -o results_nested --detection 0.00 --significant 0.00'.split(' '))

## polish muller plot using ggplot/ggmuller in R
subprocess.run('Rscript plot_polished_muller_p15.r')
# > generates figure (pdf) in working directory

# leave working subject-specific directory
os.chdir('..')

##############################
###### subject_staphAD_16
##############################

subj_id = "StaphAureus_fmk_16"
visit_time_pat = visit_time.loc[visit_time['patient'] == int(subj_id.split('_')[2])][['time1','time2','time3','time4','time5']].round(0).astype(int)

## output folder
res_folder = subj_id
os.chdir(res_folder)

## get muller/lolipop input table
lod_freq = []
ctr = 1
for index, row in data.iterrows():
    if row['0subject'] == subj_id and row['pos'] != 59855:
        # remove SNP that is fixed at t0 [5:59855]
        dd = {}
        dd['Population'] = row['0subject']
        dd['Trajectory'] = ctr
        dd['Chromosome'] = row['chr']
        dd['Position'] = row['pos']
        print(row['pos'])
        dd['Class'] = 'SNP'
        dd['Mutation'] = str(row['chr'])+"_"+str(row['pos'])+"_"+str(row['muts'])
        for j,i in enumerate(range(1,6)): # loop through visit freq data
           if row['0visit'+str(i)+'_counts'] != 0: # only record visit when samples avail (V1 always has samples!); 
                dd[str(visit_time_pat.values[0,j])] = row['0visit'+str(i)+'_freqs']
        ctr = ctr + 1
        lod_freq.append(dd)
muller_tbl = pd.DataFrame(lod_freq)
muller_tbl.to_csv("muller_tbl.tsv", sep="\t",index=False)



## run lolipop in conda base
# fix genotypes and ancestry if needed and append to input files
subprocess.run('rm -r results'.split(' '))
subprocess.run('lolipop lineage --input muller_tbl.tsv -o results --similarity-cutoff 0.001 --known-genotypes known_genotypes.csv --known-ancestry known_ancestry.tsv'.split(' ')) # sim cutoff required to correctly separate genotypes


## smooth out frequency rise by exponential function
# turn 0,1,2,3,9 > 0,100,200,300,900
# every rise in freq preceded by 10 step exponential cline (drop)

freq_arr = np.array(muller_tbl.iloc[:,6:12])
arr = freq_arr
hdr = muller_tbl.columns[6:12].astype('int') * 100


## read in trajectory > genotype array and sort by trajectory numeric id...ie. sort as mueller_tbl [ony available after first run of lolipop]
with open('results/supplementary-files/muller_tbl.genotypemembers.json') as f:
  genotype_traj_dc = json.load(f)

traj_genotype = np.array([])
for gt in genotype_traj_dc:
    for t in genotype_traj_dc[gt]:
        traj_genotype = np.append(traj_genotype,[int(t),gt])
traj_genotype = traj_genotype.reshape(int(traj_genotype.size/2),2)
traj_genotype=traj_genotype[ np.argsort(traj_genotype[:,0].astype('int')) ,:] # sort based on traj value, which turns array in same order as freq_tbl/muller_tbl


## read in edges [ony available after first run of lolipop]
edges = np.loadtxt('results/tables/muller_tbl.edges.tsv',skiprows=1,delimiter='\t',dtype=object)


## smooth freq shifts and consider nest-level for lolipop                
# extended_data_df = smooth_nested_muller_freq_table_s16(arr,hdr,edges,traj_genotype)
extended_data_df = smooth_nested_muller_freq_table(arr,hdr,edges,traj_genotype)
extended_data_df[extended_data_df > 1] = 1

extended_muller_tbl = pd.concat([muller_tbl.iloc[:,0:6], extended_data_df], axis=1, sort=False)
extended_muller_tbl.to_csv("muller_extended_tbl.tsv", sep="\t",index=False)

##run lolipop from cl
subprocess.run('rm -r results_nested'.split(' '))
subprocess.run('lolipop lineage --input muller_extended_tbl.tsv -o results_nested --detection 0.00 --significant 0.00 --similarity-cutoff 0.001 --known-genotypes known_genotypes.csv --known-ancestry known_ancestry.tsv'.split(' '))

## polish muller plot using ggplot/ggmuller in R
subprocess.run('Rscript plot_polished_muller_p16.r')
# > generates figure (pdf) in working directory

# leave working subject-specific directory
os.chdir('..')

##############################
###### subject_staphAD_26
##############################


subj_id = "StaphAureus_fmk_26"
visit_time_pat = visit_time.loc[visit_time['patient'] == int(subj_id.split('_')[2])][['time1','time2','time3','time4','time5']].round(0).astype(int)

## output folder
res_folder = subj_id
os.chdir(res_folder)


## get muller/lolipop input table
lod_freq = []
ctr = 1
for index, row in data.iterrows():
    if row['0subject'] == subj_id:
        dd = {}
        dd['Population'] = row['0subject']
        dd['Trajectory'] = ctr
        dd['Chromosome'] = row['chr']
        dd['Position'] = row['pos']
        print(row['pos'])
        dd['Class'] = 'SNP'
        dd['Mutation'] = str(row['chr'])+"_"+str(row['pos'])+"_"+str(row['muts'])
        for j,i in enumerate(range(1,6)): # loop through visit freq data
           if row['0visit'+str(i)+'_counts'] != 0: # only record visit when samples avail (V1 always has samples!); 
                dd[str(visit_time_pat.values[0,j])] = row['0visit'+str(i)+'_freqs']
        ctr = ctr + 1
        lod_freq.append(dd)
muller_tbl = pd.DataFrame(lod_freq)
muller_tbl['2'][[16]]=1 # 5_132263_['G9A'] turn freq from 0.95 to 1 at last visit. NEcessary for smoothing func, which otherwise will try to add a 6th nesting level
muller_tbl.to_csv("muller_tbl.tsv", sep="\t",index=False)



## run lolipop in conda base
# fix genotypes and ancestry if needed and append to input files
subprocess.run('rm -r results'.split(' '))
subprocess.run('lolipop lineage --input muller_tbl.tsv -o results --known-genotypes known_genotypes.csv'.split(' ')) # sim cutoff required to correctly separate genotypes

## smooth out frequency rise by exponential function
# turn 0,1,2,3,9 > 0,100,200,300,900
# every rise in freq preceded by 10 step exponential cline (drop)

freq_arr = np.array(muller_tbl.iloc[:,6:12])
arr = freq_arr
hdr = muller_tbl.columns[6:12].astype('int') * 100


## read in trajectory > genotype array and sort by trajectory numeric id...ie. sort as mueller_tbl [ony available after first run of lolipop]
with open('results/supplementary-files/muller_tbl.genotypemembers.json') as f:
  genotype_traj_dc = json.load(f)

traj_genotype = np.array([])
for gt in genotype_traj_dc:
    for t in genotype_traj_dc[gt]:
        traj_genotype = np.append(traj_genotype,[int(t),gt])
traj_genotype = traj_genotype.reshape(int(traj_genotype.size/2),2)
traj_genotype=traj_genotype[ np.argsort(traj_genotype[:,0].astype('int')) ,:] # sort based on traj value, which turns array in same order as freq_tbl/muller_tbl


## read in edges [ony available after first run of lolipop]
edges = np.loadtxt('results/tables/muller_tbl.edges.tsv',skiprows=1,delimiter='\t',dtype=object)

## smooth freq shifts and consider nest-level for lolipop                
extended_data_df = smooth_nested_muller_freq_table(arr,hdr,edges,traj_genotype)
extended_data_df[extended_data_df > 1] = 1

extended_muller_tbl = pd.concat([muller_tbl.iloc[:,0:6], extended_data_df], axis=1, sort=False)
extended_muller_tbl.to_csv("muller_extended_tbl.tsv", sep="\t",index=False)

##run lolipop from cl
subprocess.run('rm -r results_nested'.split(' '))
subprocess.run('lolipop lineage --input muller_extended_tbl.tsv -o results_nested --detection 0.00 --significant 0.00 --known-genotypes known_genotypes.csv'.split(' ')) # sim cutoff required to correctly separate genotypes

## polish muller plot using ggplot/ggmuller in R
subprocess.run('Rscript plot_polished_muller_p26.r')
# > generates figure (pdf) in working directory

# leave working subject-specific directory
os.chdir('..')

