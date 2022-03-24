#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 12:21:30 2021

@author: felix key
"""

### ***** Code for all main and SOM figures/tables ***** ###
# Input data made during previous analysis (within/across host, public data)
# Otherwise data provided (here or in metadata folder)


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
from pylab import * #what is that

from matplotlib import rc
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.font_manager import FontProperties


##########################################################################################
##########################################################################################
######### Timing dot plot, isolate quantification, and Scorad severity merged
######### Figure 1B - right panel
##########################################################################################
##########################################################################################
# dots aligned in real time
# dots show scorad severity
# dot size follows bins of # isolates

from datetime import datetime

coll_date = pd.read_table('../metadata/patient_visit_date.txt')
days_per_month = 30 # use 30d/m to convert time in days to months (entire analysis based on mnonths not days)
date_lod = []
for p in np.unique(coll_date['patient']):
    date_dc = {}
    date_dc['patient'] = p
    coll_date_p = coll_date[ (coll_date['patient'] == p) ]
    for v in [1,2,3,4,5]:
        vdate = coll_date_p[ (coll_date_p['visit'] == v) ]['date'].to_string(index=False).strip()
        if vdate == '-':
            date_dc['time'+str(v)] = np.nan
            date_dc['time'+str(v)] = np.nan
            continue
        vdate = datetime.strptime("20"+vdate.split('.')[2] + "-" + vdate.split('.')[1] + "-" + vdate.split('.')[0], '%Y-%m-%d')
        date_dc['date'+str(v)] = vdate
        if v == 1:
            date_dc['time'+str(v)] = 0
        else:
            date_dc['time'+str(v)] = (vdate - date_dc['date1']).days/days_per_month
    date_lod.append(date_dc)
pat_visit_timing = pd.DataFrame(date_lod)
pat_visit_timing.to_csv('../metadata/patient_visit_date_timepast.csv',index=False)

## order df based on isolate count and define dot plot size
tree = Phylo.read("RAxML_bestTree.Figure_1C.nwk", "newick")
pat_ctr = []
pat_visit_ctr = []
for i in tree.get_terminals():
    if i.name != 'Sref':
        pat_ctr.append(i.name.split('_')[0])
        pat_visit_ctr.append('_'.join(i.name.split('_')[:2])) 
pat_ctr = np.unique(np.array(pat_ctr),return_counts=True)
pat_srt = pat_ctr[0][np.argsort(pat_ctr[1])]
pat_srt = np.array([x[1:] for x in pat_srt],dtype=int)-1 # turn sorted list of patients in sorted list of data-df indices
pat_srt = np.insert(pat_srt,len(pat_srt)-6,14) # insert extra index, which will become white space in figure
pat_visit_timing_srt = pat_visit_timing.iloc[pat_srt[::-1],:]
pat_visit_ctr = np.unique(pat_visit_ctr,return_counts=True)


## read scoradTotal data and define dot size
meta = pd.read_csv('../metadata/subject_basic_stats_scorad.csv')
col_choice = ['#969696','#ffeda0','#feb24c','#f03b20'] # grey, YlOrRd
col_ls = []
dot_size = []
zero = 200
until10 = 201 # to allow computational distinction to zero but no visual distinction ;)
until25 = 350 # same as 50
until50 = 350
above50= 550
marker_style = []
for p in pat_visit_timing_srt['patient']:
    # scorad coloring
    for i,scr in enumerate(['1scoradtot','2scoradtot','3scoradtot','4scoradtot','5scoradtot']):
        scorad = meta[meta['Px']==p][scr].values[0]
        if isnan(scorad):
            col_ls.append(col_choice[0])
        elif scorad >= 50:
             col_ls.append(col_choice[3])
        elif scorad >= 25:
            col_ls.append(col_choice[2])
        else:
            col_ls.append(col_choice[1])
        # dot size
        if isnan(scorad):
            dot_size.append(zero)
            marker_style.append('o')
        else:
            idx_iso_ctr = np.where(pat_visit_ctr[0] == 'S'+str(p)+'_V'+str(i+1))[0]
            if pat_visit_ctr[1][idx_iso_ctr].size == 0:
                dot_size.append(zero)
                marker_style.append('o')
            elif pat_visit_ctr[1][idx_iso_ctr] <= 10:
                dot_size.append(until10)
                marker_style.append('o')
            elif pat_visit_ctr[1][idx_iso_ctr] <= 25:
                dot_size.append(until25)
                marker_style.append('o')
            elif pat_visit_ctr[1][idx_iso_ctr] <= 50:
                dot_size.append(until50)
                marker_style.append('o')
            else:
                dot_size.append(above50)
                marker_style.append('o')



y_values_pat_v = np.tile(np.arange(1,pat_visit_timing_srt.shape[0]+1)[::-1],5).reshape(5,len(pat_visit_timing_srt['patient'])).transpose()
y_values_pat_v_nonzeroIso_bool = (np.array(dot_size).reshape(len(pat_visit_timing_srt['patient']),5) != zero)
x_values_timePoints = pat_visit_timing_srt[['time1','time2','time3','time4','time5']]
x_values_timePoints_zeroOnly = np.copy(x_values_timePoints)
x_values_timePoints_zeroOnly[y_values_pat_v_nonzeroIso_bool] = np.nan

fig = plt.figure(figsize=(6,14))
plt.rc('font', family='Helvetica');plt.rcParams['pdf.fonttype'] = 42 # change font to helvetica
ax = fig.add_subplot()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

[plt.axhline(y=i, linestyle='--', color='grey', zorder=0) for i in np.arange(1,27)]
[plt.axvline(x=i, linestyle='--', color='grey', zorder=0) for i in [0,1,2,3,9]]
# plt.yaxis.grid(True,zorder=50)
# plt.vlines(x=[0,1,2,3,9],ymin=1,ymax=28,color='grey',linestyles='--')
plt.scatter(x_values_timePoints,y_values_pat_v,c=col_ls, s=dot_size,edgecolors='black')
plt.scatter(x_values_timePoints_zeroOnly,y_values_pat_v,c='black', s=zero-90,marker='x') # add x to zero samples

#np.tile(pat_visit_timing_srt['patient'],5).reshape(5,len(pat_visit_timing_srt['patient'])).transpose(),c=col_ls, s=100)
plt.xlabel("Time (months)",fontsize=14)
plt.ylabel("Patients",fontsize=14)
plt.xticks(ticks=[0,1,2,3,9,12], labels=[0,1,2,3,9,12],fontsize=18)
plt.yticks(ticks=np.arange(0,pat_visit_timing_srt.shape[0])+1,labels=pat_visit_timing_srt['patient'][::-1],fontsize=14)
plt.ylim(0,pat_visit_timing_srt.shape[0]+1)
plt.savefig('patient_visit_timing_scorad_isoCntSize.png',dpi=300,transparent=False)


## dot size legend
fig = plt.figure(figsize=(2,2))
ax = fig.add_subplot()
plt.box(False)
ax.get_yaxis().set_ticks([])
ax.get_xaxis().set_ticks([])
plt.scatter([0,0,0,0],[0,0.4,0.8,1.2],c='white', s=[200,200,450,750][::-1],edgecolors='black')
plt.savefig('patient_visit_timing_scorad_isoCntSize_legend.png',dpi=300,transparent=True)


##########################################################################################
##########################################################################################
######### Paired-stacked barplot with location and clonal complex info
######### Figure 1B - right panel
##########################################################################################
##########################################################################################
from Bio import Phylo
import numpy as np

## build source table that contoains for each isolate:
# name
# CC
# CC color 
# patient
# visit
# location

## get CC info and color
cc_cmap_dict = {'CC1':'#fbb4ae','CC5':'#fed9a6', 'CC8':'#ffffcc', 'CC15':'#ccebc5', 'CC30':'#b3cde3', 'CC45':'#decbe4', 'CC97':'#e5d8bd', 'NA':'#ffffff'} # CC1, CC5, CC8, CC15, CC30, CC45, CC97, NA
iso_st_cc = []
with open('../metadata/srst2_all_res_wSubj_wCC.txt','r') as file:
    file_wo_hdr = file.readlines()[1:] # skip header line
    for line in file_wo_hdr:
        line = line.strip().split(',')
        iso_st_cc.append([line[0],line[1],line[14],cc_cmap_dict[line[14]]])
iso_st_cc = np.array(iso_st_cc)

# read tree
tree = Phylo.read("RAxML_bestTree.Figure_1C.nwk", "newick")

## loop through tree labels and assign CC and get all info required for downstream plotting
# remove all isolates that have no identified ST (NF/failed -> low cov in ST genes). Miss some 10% of data but has no effect on proportional representation
data_iso_cc = np.array([])
iso_only = []
count_fails = 0
for i in tree.get_terminals():
    i = i.name
    if i != 'Sref':
        iso_only.append(i)
        meta = i.split('_')
        iso_dat = iso_st_cc[ iso_st_cc[:,0]==meta[3] , :]
        if iso_dat.size == 0: # should not happen!; sanity check
            print(i)
            continue
        if iso_dat[0,1] != 'failed' and 'NF' not in iso_dat[0,1]:
            data_iso_cc = np.append(data_iso_cc,np.array([meta[0],meta[1],meta[2],meta[3],iso_dat[0,2],iso_dat[0,3]]))
        else:
            count_fails = count_fails + 1
        if iso_dat[0][0][:3] == '063' or iso_dat[0][0][:3] == '073':
            # hack to include S17. Allele fingerpring not in MLST db (all isolates typed NF), but closest (one mismatch) is CC30, which is also confirmed in phylogeny. >> Manually turn S17 data to CC30
            data_iso_cc = np.append(data_iso_cc,np.array([meta[0],meta[1],meta[2],meta[3],'CC30','#b3cde3']))
data_iso_cc = data_iso_cc.reshape(int(data_iso_cc.size/6) , 6)


## write csv containing isolate-id, clonal complex, cc-color
# Note: Some CC are inferred by phylogenetic association bcs data did not allow ST typing
np.savetxt("isolate_withCCinference.csv", data_iso_cc, delimiter=",",fmt='%s')
np.savetxt("isolates_all_butRef.csv", np.array(iso_only), delimiter="\n",fmt='%s')



## remove left/right info
data_iso_cc = np.char.strip(data_iso_cc, "Left-")
data_iso_cc = np.char.strip(data_iso_cc, "Right-")
data_iso_cc[data_iso_cc[:,2] == 'Nar' , 2] = 'Nare' # strip() command removes single letter, too, thus need to fix Nare again


## build matrix for stacked barplot
cc_arr_locator = {'CC1': 0, 'CC5':1, 'CC8':2, 'CC15':3, 'CC30':4, 'CC45':5,'CC97':6, 'NA':7} # order of cc in each cc_encoder fragment in arr_cc; cc_encoder has 9th element, which is NA when NO data exists
arr_loc = np.array([]) # mx containing relative frequecny of each location for each visit (5 location, 5 visits == 25 values). loc order per visit: Pop, For, Cub, Nare, None; sum per patient 5
arr_cc = np.array([]) # same as arr_loc but for each loc 9 different cc (7xCC, NA, no isolate); cc ordered: ; 8 cc * 4 loc * 5 visits * 25 (26) patients; sum per patient 5
for pat in pat_visit_timing_srt['patient']: # use this df to maintain similar patient sorting
    pat = 'S'+str(pat)
    pat_data = data_iso_cc[ data_iso_cc[:,0]==pat , : ]
    for v in ['V1', 'V2', 'V3', 'V4', 'V5']: # loop thorugh all visits. need also info for absent visits!
        # print(pat,v)
        pat_data_v = pat_data[ pat_data[:,1]==v , : ]
        if pat_data_v.size == 0: # no isolates at visit 
            arr_loc = np.append(arr_loc,np.array([0,0,0,0,1],dtype=float)) # encode freq of loc as None [Pop, For, Cub, Nare, None]
            arr_cc = np.append(arr_cc,np.array([0,0,0,0,0,0,0,0,1],dtype=float)) # encode cc. 8th val is None
            for _ in range(3):
                arr_cc = np.append(arr_cc,np.array([0,0,0,0,0,0,0,0,0],dtype=float)) # add 3 more cc_encoder as placeholders for 3/4 remaning locations
        else:
            loc_encoder = np.array([0,0,0,0,0], dtype=float)
            for j,loc in enumerate(['Popliteal-Fossa','Forearm','Cubital-Fossa','Nare']):
                # extract visit-loc freq
                pat_data_v_loc = pat_data_v[ pat_data_v[:,2] == loc , 4 ] # only store CC info. Dont need more
                freq_loc_at_visit = pat_data_v_loc.size/pat_data_v.shape[0] # relative freq of location at visit
                loc_encoder[j] = freq_loc_at_visit
                # extract cc info; for each loc I obtain 8 cc values
                cc_encoder = np.array([0,0,0,0,0,0,0,0,0], dtype=float)
                cc_at_loc_at_visit = np.unique(pat_data_v_loc,return_counts=True)
                cc_freq_at_loc_at_visit = [cc/np.sum(cc_at_loc_at_visit[1]) for cc in cc_at_loc_at_visit[1]] # frequencies of CC at Visit-loc. Sum to 1, but will be adjusted to sum up to freq_loc_at_visit; should work for multiple CC at visit-loc
                for idx,cc_freq in enumerate(cc_freq_at_loc_at_visit):
                    cc_encoder[ cc_arr_locator[ cc_at_loc_at_visit[0][idx] ] ] =cc_freq*freq_loc_at_visit # normalize CC freq by actual freq of loc at visit; only relevant Pat21 with two CC at Nare/V3
                arr_cc = np.append(arr_cc,cc_encoder)
            arr_loc = np.append(arr_loc,loc_encoder)
# reshape: rows: visit: location;  columns: patients                
arr_loc = np.transpose(arr_loc.reshape(pat_visit_timing_srt['patient'].size, 5*5))
arr_cc = np.transpose(arr_cc.reshape(26,9*5*4))



## barplot visit data 
location_col_dc = {0:"#cccccc",1:"#969696",2:"#636363",3:"#252525",4:'white'} 
cc_cmap_list = ['#fbb4ae','#fed9a6', '#ffffcc', '#ccebc5', '#b3cde3', '#decbe4', '#e5d8bd', '#ffffff','#ffffff'] # CC1, CC5, CC8, CC15, CC30, CC45, CC97, NA
# width of the bars
barWidthR1 = 0.5
barWidthR2 = 0.3
# The x position of bars
r1 = np.arange(arr_loc.shape[1])
r2 = np.arange(arr_loc.shape[1])+0.3

fig, ax = plt.subplots(figsize=(18,7))
plt.box(False)
for i in [1,2,3,4]:
    plt.hlines(i,color='black',xmin=-0.5, xmax=25.5)
[plt.axvline(x=i+0.16, linestyle='--', color='grey') for i in r1] #, zorder=1

# add cc stacked bars
col_ctr = 0
plt.bar(r1, arr_cc[0,:], width = barWidthR1, color = cc_cmap_list[0], edgecolor = 'black')
prev_bar_data = np.copy(arr_cc[0,:])
for i in range(1,arr_cc.shape[0]):
    col_ctr += 1
    if col_ctr != 8:
        plt.bar(r1, arr_cc[i,:], bottom=prev_bar_data, width = barWidthR1, color = cc_cmap_list[col_ctr] , edgecolor = 'black')
        prev_bar_data += np.copy(arr_cc[i,:])
    else: # white edgecolor for NA cases
        plt.bar(r1, arr_cc[i,:], bottom=prev_bar_data, width = barWidthR1, color = cc_cmap_list[col_ctr] , edgecolor = 'white')
        prev_bar_data += np.copy(arr_cc[i,:])
        col_ctr = -1 # reset

# add location stacked bars
col_ctr = 0
plt.bar(r2, arr_loc[0,:], width = barWidthR2, color = location_col_dc[0], edgecolor = 'black')
prev_bar_data = np.copy(arr_loc[0,:])

for i in range(1,arr_loc.shape[0]):
    col_ctr += 1
    if col_ctr != 4:
        plt.bar(r2, arr_loc[i,:], bottom=prev_bar_data, width = barWidthR2, color = location_col_dc[col_ctr] , edgecolor = 'black')
        prev_bar_data += np.copy(arr_loc[i,:])
    else: # white edgecolor for NA cases
        plt.bar(r2, arr_loc[i,:], bottom=prev_bar_data, width = barWidthR2, color = location_col_dc[col_ctr] , edgecolor = 'white')
        prev_bar_data += np.copy(arr_loc[i,:])
        col_ctr = -1 # reset
ax.get_yaxis().set_ticks([])
ax.get_xaxis().set_ticks([])
plt.savefig('CC_visit_location_paired_stacked_barplots.png',dpi=300,transparent=True)


### get total count of patients/#isolates for ordered 
tree_names = []
for i in tree.get_terminals():
    i = i.name
    tree_names.append(i.split('_')[0])
tree_names = np.array(tree_names)
for pat in pat_visit_timing_srt['patient']: # use this df to maintain similar patient sorting
    pat = "S"+str(pat)
    print(np.sum(tree_names == pat))
    


##########################################################################################
##########################################################################################
######### Major minor clade colonising patients Pie Chart alongside USA300 tree
######### Figure 1C
##########################################################################################
##########################################################################################
# Align collapsed USA300-aligned data nodes with pie charts that have info for %samples in major or minor node
# data below extracted from tree RAxML_bestTree.Figure_1C.nwk
# _a major lineage
# _b minor lineage

cid = np.array(['22','08','13','27','26_b','20_b','17','18','04_b','25_b','23','14','21_b','07','20_a','06','02','10_a','20_c','15','21_a','09','12','11','26_a','01_b','16','10_b','19','04_a','USA300','28','01_a','25_a'])
cct = np.array([62,31,14,51,29,4,22,22,3,2,36,25,11,38,42,2,1,14,1,99,26,189,170,26,151,1,206,8,54,127,1,60,16,8])

# get freq of major minor clades
cfr = []
for i,pid in enumerate(cid):
    if len(pid.split('_')) > 1:
        tag = pid.split('_')[1]
        pid = pid.split('_')[0]
        if tag == 'a' and pid != '20':
            search_tag = 'b'
        elif pid == '20':
            # 20 special case bcs 3 lineages
            if tag=='a':
                idx=[1,2]
            elif tag=='b':
                idx=[0,2]
            else:
                idx=[0,1]
            tag_to_search = np.array(['a','b','c'])[idx]
            search_pid_idx1 = np.where(cid == pid+'_'+tag_to_search[0])[0]
            search_pid_idx2 = np.where(cid == pid+'_'+tag_to_search[1])[0]
            print(cct[i],'/',(cct[i]+cct[search_pid_idx1]+cct[search_pid_idx2])[0])
            cfr.append(cct[i]/(cct[i]+cct[search_pid_idx1]+cct[search_pid_idx2])[0])
            continue
        else:
            search_tag = 'a'
        search_pid_idx = np.where(cid == pid+'_'+search_tag)[0]
        print(cct[i],'/',(cct[i]+cct[search_pid_idx])[0])
        cfr.append( cct[i]/(cct[i]+cct[search_pid_idx])[0] )
    else:
        print(cct[i])
        cfr.append(1)
    


# plot column with pie charts. major minor pie black are complementary
plt.plot([cid.size,1,cid.size]) # 1 col
for i,clade in enumerate(cid):
    clade = clade.split('_')
    counterclock_pie = False # for plotting
    if len(clade) > 1: 
        # major/minor was defined
        if clade[1] == 'b':
            counterclock_pie = True # for plotting
    if clade[0] != 'USA300':
        plt.subplot(cid.size,1,i+1) # i+1, bcs i 0-based but not plot-counter
        plt.pie([cfr[i],(1-cfr[i])],colors=['black','white'],startangle=90,counterclock=counterclock_pie)

plt.savefig('pie_charts_usa300tree.pdf',bbox_inches='tight',transparent=True)         




####################################################################################################
####################################################################################################
######################### Pairwise Differences
######################### SOM Figure
####################################################################################################
####################################################################################################

####################################
####################################
#### Intra/Inter/major-minor pairwise diff analysis
####################################
####################################
#### Histograms of intra and inter-subject pairwise distribution separately
## build 3 arrays: 
#   (i): intra-lineage pairwise distribution (incl, major vs. major and minor vs. minor)
#   (ii) inter-lineage pairwise
#   (iii) major vs. minor lineage of subject

os.chdir('/Users/u_key/Dropbox (MIT)/Lieberman Lab/Projects/staph_aureus_AD/2021_08_1700spl_USA300/analysis/')

# required data:
# isolates major lineage only made for pan-genome assembly: /Users/u_key/Dropbox (MIT)/Lieberman Lab/Projects/staph_aureus_AD/2021_08_1700spl_USA300/analysis/isolates_major_lineages_usa300phylo_210923.txt

### plotting
## funcs
def matrix2array(mymatrix):
    ''' turn matrix to array. Only export each pair once (ie. i,j but not j,i). np.append slow... '''
    # that is terrible. there has to be a smarter approach
    pairwise_diff_array = np.array([])
    for i in np.arange(mymatrix.shape[0]):
        for j in np.arange(mymatrix.shape[1]):
            if i < j:
                pairwise_diff_array = np.append(pairwise_diff_array, mymatrix[i,j])
    return pairwise_diff_array


## data (made in analysis.py)
with open('../across_patient_analysis/pairwise_data.pickle', 'rb') as handle:
    pairwise_data = pickle.load(handle)


## turn matrix to array (remove double val) and store
# pairwise_diff_array = matrix2array(pairwise_data['pairwise_diff'])
# with open('pairwise_diff_array.pickle', 'wb') as handle:
#     pickle.dump(pairwise_diff_array, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open('pairwise_diff_array.pickle', 'rb') as handle:
    pairwise_diff_array = pickle.load(handle)


## get array with isolates of major/minor lineages
treesampleNamesLong = pairwise_data['treesampleNamesLong']
all_major_lineage_isolates = np.loadtxt('../metadata/isolates_major_lineages_usa300phylo_210923.txt',delimiter=",",dtype=object)[:,1] # only major lineage isolates

# obtain all isolates from minor lineages
maj_min_info = [] # var storing all isolates of subjects with major/minor lineages. Subjects w/ single lineage not included!
for s in treesampleNamesLong:
    kit_spl = s.split('_')[3]
    if kit_spl not in all_major_lineage_isolates:
        maj_min_info.append([s,'Minor'])
maj_min_info = np.array(maj_min_info)

# add all isolates of major lineages for patients with identified minor lineage
maj_min_subject = np.unique(np.array([ i.split('_')[0] for i in maj_min_info[:,0]]))
for s in treesampleNamesLong:
    if s.split('_')[0] in pat_w_minor:
        if s not in maj_min_info[:,0]:
            maj_min_info = np.vstack([maj_min_info,[s,"Major"]])

maj_min_info[np.where(maj_min_info[:,0]=='S20_V1_Right-Popliteal-Fossa_069-RD9')[0],1] = "Minor2"

## get pairwise count matrices for different groups. 
# pairwise groups only considered if 75% of sites were called in both isolates.
min_p_usable_count = 0.75 * 60989 # 60989 is total # of SNPs in USA300 tree

pairwise_diff = pairwise_data['pairwise_diff']
subjectID = [i.split('_')[0] for i in treesampleNamesLong]


pairwise_intra_array = np.array([])
pairwise_inter_array = np.array([])
pairwise_inter_array_S4S19 = np.array([])
pairwise_majmin_array = np.array([])
closely_related_interlineage_ids = np.array([])
for i in np.arange(pairwise_diff.shape[0]):
    if pairwise_data['p_usable'][i] < min_p_usable_count:
        continue
    for j in np.arange(pairwise_diff.shape[1]):
        if pairwise_data['p_usable'][j] < min_p_usable_count:
            continue
        if i < j: # does not make a difference to restrict on usable sites (see 50pofpusable_* pdfs): and pairwise_data['p_usable'][i] > 35429 and pairwise_data['p_usable'][j] > 35429 : # at least half of p was usable (avoid artefaficially low pairwise diff due to low cov!), 70857/2 = 35429
            if subjectID[i] == subjectID[j]: # pair same subject
                if subjectID[i] in maj_min_subject: # it is a subject that has a major/minor lineage
                    if maj_min_info[ maj_min_info[:,0]==treesampleNamesLong[i],1 ] == maj_min_info[ maj_min_info[:,0]==treesampleNamesLong[j],1 ]: # it is major/major or minor/minor
                        pairwise_intra_array = np.append(pairwise_intra_array, pairwise_diff[i,j])
                    else:
                        pairwise_majmin_array = np.append(pairwise_majmin_array, pairwise_diff[i,j])
                else: # intra diversity of subjects without major/minor
                    pairwise_intra_array = np.append(pairwise_intra_array, pairwise_diff[i,j]) 
            else:
                subject_pair = subjectID[i]+subjectID[j]
                if subject_pair == 'S19S4' or subject_pair == 'S4S19':
                    pairwise_inter_array_S4S19 = np.append(pairwise_inter_array_S4S19, pairwise_diff[i,j])
                else:
                    pairwise_inter_array = np.append(pairwise_inter_array, pairwise_diff[i,j])
                if pairwise_diff[i,j] < 100:                    
                    closely_related_interlineage_ids = np.append(closely_related_interlineage_ids, (treesampleNamesLong[i],treesampleNamesLong[j]) )
                    if subjectID[i] != 'S4' and subjectID[i] != 'S4':
                        print(i,j,treesampleNamesLong[i],treesampleNamesLong[j])
                        print(pairwise_data['p_usable'][i],pairwise_data['p_usable'][j])


## plot intra patient pairwise diversity
bins=np.exp2(np.linspace(0, 15, 60))
fig,ax1 = plt.subplots(figsize=(4,3))
ax1.hist(pairwise_intra_array, bins, density=False, facecolor='#d95f02', alpha=0.75,label='Intra patient pairwise differences')
ax1.hist(pairwise_majmin_array, bins, density=False, facecolor='#d95f02', alpha=0.75,label='Major vs. minor lineage')
plt.xlabel('Pairwise Distances (SNVs)')
plt.ylabel('Counts')
plt.yscale('log', nonposy='clip')
plt.xscale('log')
plt.xlim(1,30000)
# plt.legend(loc='upper right',fontsize='small')
plt.savefig('pairwise_intra_majmin_USA300FPR3757_mapped_log.pdf',bbox_inches = "tight")
plt.close()

## plot inter patient pairwise distance
bins=np.exp2(np.linspace(0, 15, 60))
fig,ax1 = plt.subplots(figsize=(4,3))
n, bins, patches = ax1.hist(pairwise_inter_array, bins, density=False, facecolor='#7570b3', alpha=0.75,label='Inter lineage pairwise differenes')
ax1.hist(np.array([0]), bins,  density=False, facecolor='#cb181d', alpha=0.5,label='S04-S19 pairwise differences') # visually differentiate S4-S19 pairs vs all other
plt.xlabel('Pairwise Distances (SNVs)')
plt.ylabel('Counts')
plt.yscale('log', nonposy='clip')
plt.xscale('log')
plt.ylim(1, 8000000)
plt.xlim(1, 30000)
plt.savefig('pairwise_inter_USA300FPR3757_mapped_log.pdf',bbox_inches = "tight")
plt.show()
plt.close()


############
###### SOM table 1. Minor-Major lineage overview
############
d=[]
for r in maj_min_info:
    iso=r[0].split('_')
    new_string=iso[0]+'_'+iso[1]+'_'+iso[2]+'_'+r[1]
    d.append(new_string)
d=np.array(d)

res=np.unique(np.sort(d),return_counts=True)
dd=[]
for i,r in enumerate(res[0]):
    rs=r.split('_')
    dd.append( [ rs[0],rs[3],rs[1],rs[2],res[1][i] ] )
dd=np.array(dd)
np.savetxt('major_minor_overview.csv',dd,delimiter=',',fmt='%s')
    
    
    
##########################################################################################
##########################################################################################
######### Scatter & Correlation #-isolates per patient visit vs. Scorad
######### Scatter & correlation %relative abundance of aureus at a visit per site versus the number of isolates visit
######### Plot fraction of visits with more than 1 lineage
######### Extended Data Figure 2 ABC
##########################################################################################
##########################################################################################
# get isolate count from USA300 tree
# get scoradTot per V from table

from Bio import Phylo
import scipy

tree = Phylo.read("RAxML_bestTree.Figure_1C.nwk", "newick")

sub_visit = np.array([],dtype=object)
for leaf in tree.get_terminals():
    if not leaf.name == 'Sref': 
        sub=leaf.name.split('_')[0]
        visit=leaf.name.split('_')[1]
        if not sub == 'S2-control':
            sub_visit = np.append(sub_visit,sub+"_"+visit)
    

sub_visit = np.unique(sub_visit,return_counts=True)

scoradData = np.genfromtxt('../metadata/subject_basic_stats_scorad.csv',dtype=float,delimiter=",",skip_header=True)[:, [0,5,7,9,11,13]] # Px,edad,genero,cloro si/no,1scoradobj,1scoradtot,2scoradobj,2scoradtot,3scoradobj,3scoradtot,4scoradobj,4scoradtot,5scoradobj,5scoradtot,MeanObj,MeanTot

x_num_isolates = np.array([],dtype=int)
y_scorad = np.array([],dtype=float)
for i in np.arange(scoradData.shape[0]): # iteration subjects
    for j in np.arange(1,6): # iteration visits
        if not np.isnan(scoradData[i,j]): # were samples collected?
            y_scorad = np.append(y_scorad,scoradData[i,j])
            s_v = 'S'+str(int(scoradData[i,0]))+'_V'+str(j) # mirror string that represents Patient_Visit in sub_visit
            idx = np.where(sub_visit[0] == s_v)
            if idx[0].size == 0:
                x_num_isolates = np.append(x_num_isolates,0)
            else:
                x_num_isolates = np.append(x_num_isolates,sub_visit[1][idx])


slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x_num_isolates , y_scorad )
# LinregressResult(slope=0.6922562531250595, intercept=25.187395144438355, rvalue=0.6033125751090774, pvalue=6.066704838595425e-14, stderr=0.0818470101609847)
regress_dict = {}
regress_dict['slope'] = slope
regress_dict['intercept'] = intercept
regress_dict['r_sq'] = r_value**2
regress_dict['p_val'] = p_value
regress_dict['std_err'] = std_err
# No CI calculated. All wyas I found required gaussian dist. Not an assumption I necessarily belive is true. Bootstrap prob best.

plt.rc('font', family='Helvetica');plt.rcParams['pdf.fonttype'] = 42 # change font to helvetica
plt.plot( x_num_isolates,y_scorad, linestyle='none', marker='o',color='grey')
plt.plot( [np.min(x_num_isolates),np.max(x_num_isolates)],(slope*np.array([np.min(x_num_isolates),np.max(x_num_isolates)])+intercept),linewidth=3,color='black')
plt.xlabel('Number isolates per Visit',fontsize=14)
plt.ylabel('SCORAD at Visit',fontsize=14)
plt.savefig('linreg_scorad_isolatesPerVisit.png',dpi=300)
plt.show()
    
##########################################################################################
######### Scatter & Correlation #-isolates per patient visit vs. 16S %SA
######### SOM
##########################################################################################
# recycle Khadka, Key et al. 2021 
asv = pd.read_csv('16s_results_relative_freq_saureus.csv') # data from VK early '22. final analysis (published)
tree = Phylo.read("RAxML_bestTree.Figure_1C.nwk", "newick")

iso = []
for leaf in tree.get_terminals():
    iso.append(leaf.name)
iso =np.array(iso) 

# build vector 16S %SA and vector #isolates
percSA=[]
numIso=[]
for index, row in asv.iterrows():
    percSA.append(row['Rel_Staph_Aureus'])
    pat,bla,v = row['SampleName'].split('_')
    number=np.sum([s.startswith('S'+pat+'_V'+v) for s in iso])
    numIso.append(number)
    print(row['SampleName'],number,row['Rel_Staph_Aureus'])
    

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(numIso,percSA)

plt.rc('font', family='Helvetica');plt.rcParams['pdf.fonttype'] = 42 # change font to helvetica
plt.plot( numIso , percSA, linestyle='none', marker='o',color='grey')
plt.plot( [np.min(numIso),np.max(numIso)],(slope*np.array([np.min(numIso),np.max(numIso)])+intercept),linewidth=3,color='black')
plt.xlabel('Number isolates per visit',fontsize=14)
plt.ylabel('Relative abundance of S.aureus per visit',fontsize=14)
plt.savefig('linreg_16s_isolateCount_pervisit.png',dpi=300)
plt.show()




##########################################################################################
######### Plot fraction of visits with more than 1 lineage
######### Cumulative plot of visits with x isolates to fraction of visits with >1 lineage
##########################################################################################
pat_visit_timing = pd.read_csv('/Users/u_key/Documents/mit/stapAD/tables/metadata/patient_visit_date_timepast.csv')
## get counts of isolates per visit across all subjects
tree = Phylo.read("/Users/u_key/Dropbox (MIT)/Lieberman Lab/Projects/staph_aureus_AD/2021_08_1700spl_USA300/analysis/run2_partialDel005_minorAF003/RAxML_bestTree.Figure_1C.nwk", "newick")
pat_visit_ctr = []
for i in tree.get_terminals():
    if i.name != 'Sref':
        pat_visit_ctr.append('_'.join(i.name.split('_')[:2]))
pat_visit_ctr = np.unique(pat_visit_ctr,return_counts=True)

# list of visits with >1 lineage
patient_visit_with_multiple_lineages = ['S4_V1','S4_V4','S21_V3','S25_V1','S20_V1']

isolatesVisit_fractionMultLin_cumulativeVisits = []
for i in range(2,71,1):
    visits_to_consider = pat_visit_ctr[1] < i
    num_visits_with_multiple_lineages = np.sum([True for j in patient_visit_with_multiple_lineages if j in pat_visit_ctr[0][visits_to_consider]])
    isolatesVisit_fractionMultLin_cumulativeVisits.append([i,num_visits_with_multiple_lineages/np.sum(visits_to_consider),np.sum(visits_to_consider)/len(visits_to_consider)])
isolatesVisit_fractionMultLin_cumulativeVisits = np.array(isolatesVisit_fractionMultLin_cumulativeVisits)


#define colors to use
col1 = 'black'# 'steelblue'
col2 = 'grey'

plt.rc('font', family='Helvetica');plt.rcParams['pdf.fonttype'] = 42 # change font to helvetica
fig,ax = plt.subplots()
ax.plot(isolatesVisit_fractionMultLin_cumulativeVisits[:,0],isolatesVisit_fractionMultLin_cumulativeVisits[:,1], color=col1,linewidth=2)
ax.set_xlabel('Number isolates per visit < x\n(cumulative)', fontsize=14)
ax.set_ylabel('Fraction of visits with\nmultiple lineages', color=col1, fontsize=16)
plt.tight_layout(pad=1.1, w_pad=0.5, h_pad=0.5)
fig.savefig('/Users/u_key/Documents/mit/stapAD/pdf/collectors_curve/fraction_visits_multiple_lineages.png',dpi=300)#,transparent=True)





##########################################################################################
##########################################################################################
######### Molecular clock 
######### Figure 2B
##########################################################################################
##########################################################################################
# mutations per isolate corrected for number of genome-wide base calls
# Incl. correct patient-specific visiting times
mol_clock_data_dc = pickle.load(open( '../within_patient_analysis/molecular_clock_data.pk1','rb'))

scale_factor_per_base = 1000000 # scales all rates by x number of bases
tinv = lambda p, df: abs(stats.t.ppf(p/2, df)) # func2turn SE2ConfInterval; source https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html
confint = 0.95

linreg_res = {}
for s in ['StaphAureus_fmk_4','StaphAureus_fmk_9','StaphAureus_fmk_12','StaphAureus_fmk_15','StaphAureus_fmk_16','StaphAureus_fmk_26']:
    x = mol_clock_data_dc[s][:,0]
    y = np.array(mol_clock_data_dc[s][:,1]/mol_clock_data_dc[s][:,2]) * scale_factor_per_base
    slope,intercept,rval,pval,stderr = stats.linregress(x=x,y=y)
    ts = tinv(1-confint, len(x)-2) 
    linreg_res[s] = [slope,intercept,(stderr*ts),(stderr*ts) ] # calc CI for slope 
    print(slope,stderr,(stderr*ts))
    #linreg_res[s] = [slope,intercept,slope-(stderr*ts),slope+(stderr*ts) ] # calc CI for slope 
    
### linreg plot molecular rate as barplot
myfontsize=18
rate_per_year = [linreg_res[x][0]*12 for x in linreg_res.keys()]
ci_yerrbar = np.vstack((np.array([linreg_res[x][2]*12 for x in linreg_res.keys()]),np.array([linreg_res[x][3]*12 for x in linreg_res.keys()])))
higher_lit_clock_rate=3.3e-6*2800000/2.8 # based on published rates for S.aureus
lower_lit_clock_rate=1.2e-6*2800000/2.8
with plt.style.context('default'):
    plt.rc('font', family='Helvetica');plt.rcParams['pdf.fonttype'] = 42 # change font to helvetica
    fig, ax = plt.subplots(figsize=(6, 4))
    plt.errorbar(np.arange(1,7),rate_per_year,color='black',yerr=ci_yerrbar,fmt='o')
    plt.axhspan(higher_lit_clock_rate,lower_lit_clock_rate, color='gainsboro', alpha=0.75, lw=0)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xticks(np.arange(1,7), [x.split('_')[2] for x in linreg_res.keys()],fontsize=myfontsize)
    plt.yticks(fontsize=myfontsize)
    plt.xlabel('Patient',fontsize=myfontsize)
    plt.ylabel('Mutation rate per 1Mb/year',fontsize=myfontsize+2)
    plt.savefig('barplot_linreg_molecular_clock.png',dpi=300)
    plt.show()


   
####################################################################################################
####################################################################################################
######################### tMRCA - @ V1 all subjects w/ > 10 isolates at V1 -- scatter plot
######################### Figure 2K
####################################################################################################
####################################################################################################
import matplotlib.patches as mpatches # for legend plotting control

## distance for each isolate obtained (saved) in analysis.py
data = pickle.load(open( '../within_patient_analysis/molecular_clock_data.pk1','rb'))

## load age of subjects
patient_age_dc = {}
with open('../metadata/age_patients_staphAD.txt') as file:
    for line in file:
        line = line.strip().split('\t')
        patient_age_dc['StaphAureus_fmk_'+line[0]] = int(line[1])

# cols
my_cols = {'StaphAureus_fmk_4':'#de2d26','StaphAureus_fmk_9':'#d95f0e','StaphAureus_fmk_12':'#2c7fb8','StaphAureus_fmk_15':'#ffd92f','StaphAureus_fmk_16':'#c51b8a','StaphAureus_fmk_26':'#31a354'}
scale_factor_per_base = 1000000
linreg_res_median = np.median(np.array([linreg_res[s][0] for s in linreg_res.keys()]))

## get data for plotting in 4 lists:
age_x = [] #...age in y
mean_scaled_dist_per_visit_mutcorr = np.array([]) # mean genetic distance based on scaled num mut
area_s = [] # define area due to subject (make highly colonized subjects bigger)
col_v = [] #...color based on visit
subj_ls = [] # subject list just for dbl check
grey_scale_idx = 0
myfontsize=18
for s in data.keys():
    if not s == 'StaphAureus_fmk_2-control':
        s_data = data[s]
        ## get first visit with 10 or more isolates
        v_isolate_counts = np.unique(s_data[:,0],return_counts=True) # also sorts!
        for j,num in enumerate(v_isolate_counts[1]):
            if num >= 10:
                v = v_isolate_counts[0][j]
                # print(s,v)
                break
        # for v in 0: #np.unique(s_data[:,0]):
        v_s_data = s_data[ np.where(s_data[:,0] == v) ]
        if v_s_data.shape[0] >= 10:
            age_x.append(patient_age_dc[s])
            subj_ls.append(s)
            if s in my_cols.keys():
                col_v.append(my_cols[s])
                area_s.append(100)
                mean_scaled_dist_per_visit_mutcorr = np.append(mean_scaled_dist_per_visit_mutcorr,((np.mean(v_s_data[:,1]/v_s_data[:,2])* scale_factor_per_base)/linreg_res[s][0])) # divide scaled de novo mut by mut rate per month
                print(s,mean_scaled_dist_per_visit_mutcorr[-1])
            else:
                # col_v.append(grey_scale[grey_scale_idx])
                col_v.append('black')
                area_s.append(50)
                mean_scaled_dist_per_visit_mutcorr = np.append(mean_scaled_dist_per_visit_mutcorr,((np.mean(v_s_data[:,1]/v_s_data[:,2])* scale_factor_per_base)/linreg_res_median)) # divide scaled de novo mut by mut rate per month
                print(s,mean_scaled_dist_per_visit_mutcorr[-1])
                grey_scale_idx += 1
myfontsize=14
plt.style.use('default')
plt.rc('font', family='Helvetica');plt.rcParams['pdf.fonttype'] = 42 # change font to helvetica
plt.figure(figsize=(6,3))
plt.scatter(age_x, mean_scaled_dist_per_visit_mutcorr, s=area_s, c=col_v, alpha=1)
plt.xlabel("Age of patient (years)",fontsize=myfontsize)
plt.ylabel("tMRCA (years)",fontsize=myfontsize)
plt.ylim(0,24)
plt.yticks(ticks=[0,12,24], labels=[0,1,2],fontsize=myfontsize)
plt.xticks(ticks=[6,8,10,12,14], labels=[6,8,10,12,14],fontsize=myfontsize)
handles_leg = {}
for i in range(len(age_x)):
    id_short = subj_ls[i].split('_')[-1]
    if len(id_short) == 1:
        id_short = "Patient 0"+id_short
    else:
        id_short = "Patient "+id_short
    handles_leg[id_short] = mpatches.Patch(color=col_v[i], label=id_short)
handles_leg['other'] = mpatches.Patch(color='black', label='other Patients')
plt.legend(handles=[handles_leg[j] for j in ['Patient 04','Patient 09','Patient 12','Patient 15','Patient 16','Patient 26','other']],loc='upper left',ncol=1,fontsize=myfontsize-4,frameon=False)
plt.savefig('gendistance_ageS_min10isol_tmrca_scaled_mutRateCorrinMonths_first10isoVisit.png',dpi=300)
plt.show()


    
##########################################################################################
##########################################################################################
######### Pie chart / donut plot capsule +/- subject colonization
######### Figure 2
##########################################################################################
##########################################################################################
# patients colonized by cap- lineage (based on Jean Lee serotyping results)
# 1 (major and minor)
# 10 (major)
# 12
# 14
# 16
# 25 (major)
# 28

## > 7/25

## three categories
# Acquired sweeping mutation in capsule locus -- S9/S12
# Ancestral stop mutation in capsule locus -- lineages: 7
# Lineages with maintained full capsule locus --- 30-8 = 22

import matplotlib.pyplot as plt
plt.rc('font', family='Helvetica');plt.rcParams['pdf.fonttype'] = 42 # change font to helvetica


# Pie chart, where the slices will be ordered and plotted counter-clockwise:
labels = 'Acquired sweeping mutation', 'Ancestral capsule loss','Maintained complete locus'
sizes = [2,7,22]

fig1, ax1 = plt.subplots()
wedges, texts = ax1.pie(sizes, labels=None, autopct=None,
        shadow=False, startangle=90, colors=['dimgray','silver','#f0f0f0']) # colors=['#3182bd','#9ecae1','#f0f0f0']
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
ax1.legend(wedges, labels,title='Capsule status', loc="center left")
fig1.savefig("pie_sweep_loss_func_v2.svg",transparent=True)
plt.show()




##########################################################################################
##########################################################################################
######### Collectors curve SNPVs per patient per visit
######### Extended Data Figure X 
##########################################################################################
##########################################################################################
import random

# approach
# load snp_table.csv
# loop through visits (min # iso per visit -- 10 (like tMRCA!))
# get mean value (95CI) of SNVs per cumulative isolate (100 permutations) vs Ref (pangenome)
# plot (sgl plot per visit?)

num_permutations = 100


## simulation #SNPs identified
# 100 times randomly choose x isolates (0<x<n) and assess mean / SD numbers SNVs identified 
# identified versus nt_anc, based on single-iolate nt_anc (if . > discarded)
# mean for 100sims, SD for 100 means

# TMRCA
res_snps = {}
for p in range(1,29):
    p = str(p)
    if not os.path.exists('../within_patient_analysis/subject_'+p+'/snp_table.csv'):
        continue # no data
    print(p)
    d=pd.read_csv('../within_patient_analysis/subject_'+p+'/snp_table.csv')
    anc_true_bool = (np.array(d[['nt_anc']]) != '?')[:,0] # bool to remove sites w/o inferred ancestral allele    
    nt_ref_arr = np.array(d[['nt_ref']])[anc_true_bool,:]
    nt_alt_arr = np.array(d[['nt_alt']])[anc_true_bool,:]
    nt_anc_arr = np.array(d[['nt_anc']])[anc_true_bool,:]
    nt_der_arr = np.array(d[['nt_alt']])[anc_true_bool,:]
    nt_der_arr[ nt_ref_arr != nt_anc_arr ] = nt_ref_arr[nt_ref_arr != nt_anc_arr]    
    for v in ['1','2','3','4','5']:
        dv = d.filter(regex=('_V'+v+'_')) # df visit
        arr = np.array([])
        if dv.empty or len(dv.columns) < 10:
            res_snps[p+'_'+v] = arr
        else:
            num_isolates = dv.shape[1]+1 # +1 for 0-based range for loop and to capture full data (all isolates) 
            for i in range(1,num_isolates): # loop through sets of isolates
                for _ in range(num_permutations): # loop permutations
                    selection_isolates = random.choices(list(dv.columns),k=i)  # random selection of i isolates; w/o replacement
                    dvi = np.array(dv[ selection_isolates ])[anc_true_bool,:]  # dataframe for i isolates of all polarized sites 
                    nt_anc_isolates = np.any((dvi == nt_anc_arr),axis=1) # any pos where any isolate has ref
                    nt_der_isolates = np.any((dvi == nt_der_arr),axis=1) # any pos where any isolate has alt
                    nt_polymorphic_sites = np.all((nt_anc_isolates,nt_der_isolates),axis=0) # pos where alt and ref present, ie. polymorphic among selected isolates
                    num_variable_snv_found_in_subset = np.sum(nt_polymorphic_sites)
                    arr = np.append(arr,num_variable_snv_found_in_subset) # store # diffs
            res_snps[p+'_'+v] = arr.reshape((dvi.shape[1],num_permutations)) 


# DMRCA
res_mean_dmrca = {}
res_sd_dmrca = {}
for p in range(1,29):
    p = str(p)
    if not os.path.exists('../within_patient_analysis/subject_'+p+'/snp_table.csv'):
        continue # no data
    print(p)
    d=pd.read_csv('../within_patient_analysis/subject_'+p+'/snp_table.csv')
    anc_true_bool = (np.array(d[['nt_anc']]) != '?')[:,0] # bool to remove sites w/o inferred ancestral allele    
    nt_ref_arr = np.array(d[['nt_ref']])[anc_true_bool,:]
    nt_alt_arr = np.array(d[['nt_alt']])[anc_true_bool,:]
    nt_anc_arr = np.array(d[['nt_anc']])[anc_true_bool,:]
    nt_der_arr = np.array(d[['nt_alt']])[anc_true_bool,:]
    nt_der_arr[ nt_ref_arr != nt_anc_arr ] = nt_ref_arr[nt_ref_arr != nt_anc_arr]    
    for v in ['1','2','3','4','5']:
        arr_mean_dmrca = np.array([])
        arr_sd_mean_dmrca = np.array([])
        dv = d.filter(regex=('_V'+v+'_')) # df visit
        if dv.empty or len(dv.columns) < 10:
            res_mean_dmrca[p+'_'+v] = arr_mean_dmrca
            res_sd_dmrca[p+'_'+v] = arr_sd_mean_dmrca
        else:
            num_isolates = dv.shape[1]+1 # +1 for 0-based range for loop and to capture full data (all isolates) 
            for i in range(1,num_isolates): # loop through sets of isolates
                isolate_iteration_dmrca_mean = np.array([])
                for _ in range(num_permutations): # loop permutations
                    selection_isolates = random.choices(list(dv.columns),k=i)  # random selection of i isolates; w/o replacement
                    dvi = np.array(dv[ selection_isolates ])[anc_true_bool,:]  # dataframe for i isolates of all polarized sites 
                    nt_anc_isolates = np.any((dvi == nt_anc_arr),axis=1) # any pos where any isolate has ref
                    nt_der_isolates = np.any((dvi == nt_der_arr),axis=1) # any pos where any isolate has alt
                    nt_polymorphic_sites = np.all((nt_anc_isolates,nt_der_isolates),axis=0) # pos where alt and ref present, ie. polymorphic among selected isolate
                    dvi_var = dvi[nt_polymorphic_sites,] # data only for variable positions in subset of isolates
                    dmrca_mean = np.mean(np.sum(dvi_var == nt_der_arr[nt_polymorphic_sites],axis=0)) # collect mean dMRCA for each iteration/isolate
                    isolate_iteration_dmrca_mean = np.append(isolate_iteration_dmrca_mean,dmrca_mean) # store mean dMRCAs
                arr_mean_dmrca = np.append(arr_mean_dmrca,np.mean(isolate_iteration_dmrca_mean)) # single value for each iteration
                arr_sd_mean_dmrca = np.append(arr_sd_mean_dmrca,np.std(isolate_iteration_dmrca_mean)) # single value for each iteration based on each mean obs 
            res_mean_dmrca[p+'_'+v] = arr_mean_dmrca
            res_sd_dmrca[p+'_'+v] = arr_sd_mean_dmrca

## plotting dMRCA
visit_col_dc = {"V1":"#BCBDDC","V2":"#9E9AC8","V3":"#807DBA","V4":"#6A51A3","V5":"#3F007D"}
myfontsize = 20

fig = plt.figure(figsize=(25,20))
plt.rc('font', family='Helvetica');plt.rcParams['pdf.fonttype'] = 42 # change font to helvetica
ctr_subplots = 1
visit_plot = {}
for p in np.unique([int(k.split('_')[0]) for k in res_mean_dmrca.keys()]):
    p = str(p)
    if p in ['10','6']:
        continue # skip patients with no visit >10 isolates
    ax = fig.add_subplot(5,4,ctr_subplots)
    for v in [1,2,3,4,5]:
        res_id = p+'_'+str(v)
        if np.any(res_mean_dmrca[res_id]):
            visit_plot[str(v)] = ax.errorbar(range(1,res_mean_dmrca[res_id].shape[0]+1), res_mean_dmrca[res_id], res_sd_dmrca[res_id]/2, linestyle='None', marker='o',color=visit_col_dc['V'+str(v)])
    # ax.set_xlabel('Number of isolates',fontsize=14)
    if ctr_subplots in [1,5,9,13,17]:
        ax.set_ylabel('dMRCA',fontsize=myfontsize)
    if ctr_subplots > 16:
        ax.set_xlabel('Number of isolates',fontsize=myfontsize)
    ax.set_title('Patient '+p,loc='left',fontsize=myfontsize)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=myfontsize)
    ctr_subplots += 1
fig.set_tight_layout({"pad": .5})
fig.savefig('dMRCAonly_allPw10iso_'+str(num_permutations)+'sim_wReplacement_v1.pdf',dpi=300)        




### plotting SNPs
fig = plt.figure(figsize=(25,20))
plt.rc('font', family='Helvetica');plt.rcParams['pdf.fonttype'] = 42 # change font to helvetica
ctr_subplots = 1
visit_plot = {}
for p in np.unique([int(k.split('_')[0]) for k in res_snps.keys()]):
    p = str(p)
    if p in ['10','6']:
        continue # skip patients with no visit >10 isolates
    ax = fig.add_subplot(5,4,ctr_subplots)
    
    for v in [1,2,3,4,5]:
        res_id = p+'_'+str(v)
        if np.any(res_snps[res_id]):
            visit_plot[str(v)] = ax.errorbar(range(1,res_snps[res_id].shape[0]+1), np.mean(res_snps[res_id],axis=1), np.std(res_snps[res_id],axis=1)/2, linestyle='None', marker='o',color=visit_col_dc['V'+str(v)])
    if ctr_subplots in [1,5,9,13,17]:
        ax.set_ylabel('SNPs identified',fontsize=myfontsize)
    if ctr_subplots > 16:
        ax.set_xlabel('Number of isolates',fontsize=myfontsize)
    ax.set_title('Patient '+p,loc='left',fontsize=myfontsize)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=myfontsize)
    ctr_subplots += 1
    
    
    
# fig.legend([visit_plot['1'],visit_plot['2'],visit_plot['3'],visit_plot['4'],visit_plot['5']],['Visit 1','Visit 2','Visit 3','Visit 4','Visit 5'],loc='lower center')
fig.set_tight_layout({"pad": .5})
fig.savefig('SNPonly_allPw10iso_'+str(num_permutations)+'sim_wReplacement_v1.pdf',dpi=300)        



##########################################################################################
##########################################################################################
######### Mutation Spectrum per patient 
######### Extended Data Figure X
##########################################################################################
##########################################################################################
# Show top 6 Patients individually. All other patients with at least 10 on-person SNVs are combined.
# SNVs polarized by ancestral nucleotide (patient specific outgroup)
nt = np.array(['A','T','C','G'])
allmuts = np.array(['AT','TA','AC','TG','AG','TC','GT','CA','GC','CG','CT','GA']); # 'AT' denotes A to T
muts = pd.read_csv('../within_patient_analysis/allGoodpos_allSubjects.csv') 
## get all mut counts in ordered (allmuts) fashion
mut_ctr = {}
for s in  np.unique(muts['subject']):
    # nts = np.array(muts.loc[ muts['subject'] == s , 'nts']) # this is ref/alt, which does loose the polarization which is important
    muts_sub = muts.loc[ muts['subject'] == s , :]
    nts = []
    for i,a in enumerate(muts_sub['nt_anc']):
        if a in nt:
            if muts_sub['nt_alt'].iloc[i] in nt and muts_sub['nt_alt'].iloc[i] != a:
                # print(a,muts_sub['der'].iloc[i])
                nts.append(a+muts_sub['nt_alt'].iloc[i])
            elif muts_sub['nt_ref'].iloc[i] in nt and muts_sub['nt_alt'].iloc[i] == a:
                # print(i,a,muts_sub['ref'].iloc[i])
                nts.append(a+muts_sub['nt_ref'].iloc[i])
    nts = np.array(nts)            
    if nts.size >= 10: # at least 10 Mutations
        mut_ctr[s] = np.zeros(len(allmuts))
        nts_ctr = np.unique(nts,return_counts=True)  
        for i,n in enumerate(nts_ctr[0]):
            # print(i)
            if len(n) == 2: # exclude triallelic
                allmuts_idx = np.where(allmuts==n)[0][0]
                mut_ctr[s][allmuts_idx] = nts_ctr[1][i]

## get vars for grouped barplots
# top5, and all other merged
s='StaphAureus_fmk_4'; s4 = [ np.sum( mut_ctr[s][i:(i+2)])/np.sum(mut_ctr[s])  for i in np.arange(0,12,2)]
s='StaphAureus_fmk_9'; s9 = [ np.sum( mut_ctr[s][i:(i+2)])/np.sum(mut_ctr[s])  for i in np.arange(0,12,2)]
s='StaphAureus_fmk_12'; s12 = [ np.sum( mut_ctr[s][i:(i+2)])/np.sum(mut_ctr[s])  for i in np.arange(0,12,2)]
s='StaphAureus_fmk_15'; s15 = [ np.sum( mut_ctr[s][i:(i+2)])/np.sum(mut_ctr[s])  for i in np.arange(0,12,2)]
s='StaphAureus_fmk_16'; s16 = [ np.sum( mut_ctr[s][i:(i+2)])/np.sum(mut_ctr[s])  for i in np.arange(0,12,2)]
s='StaphAureus_fmk_26'; s26 = [ np.sum( mut_ctr[s][i:(i+2)])/np.sum(mut_ctr[s])  for i in np.arange(0,12,2)]
# combine all other subjects
sa=['StaphAureus_fmk_11','StaphAureus_fmk_14','StaphAureus_fmk_22','StaphAureus_fmk_23','StaphAureus_fmk_27','StaphAureus_fmk_7','StaphAureus_fmk_8','StaphAureus_fmk_20','StaphAureus_fmk_19']; 
for s in sa:
    if s == 'StaphAureus_fmk_11': # first iteration
        t = np.array([ np.sum( mut_ctr[s][i:(i+2)])/np.sum(mut_ctr[s])  for i in np.arange(0,12,2)])
    else:
        t = np.concatenate((t,np.array([ np.sum( mut_ctr[s][i:(i+2)])/np.sum(mut_ctr[s])  for i in np.arange(0,12,2)])))
so = t.reshape(9,6)
so = np.mean(so,axis=0)

mut_res = np.array(s4+s9+s12+s15+s16+s26+list(so)).reshape(7,6)

# Set position of bar on X axis
barWidth = 0.1
r1 = np.arange(7)
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]
r4 = [x + barWidth for x in r3]
r5 = [x + barWidth for x in r4]
r6 = [x + barWidth for x in r5]
r7 = [x + barWidth for x in r6]
 
my_cols = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f']
barWidth=0.1
myfonts=20
# Make the plot
plt.figure(figsize=(20,5))
plt.bar(r1, mut_res[:,0], color=my_cols[0], width=barWidth, edgecolor='white', label='A->T,T->A')
plt.bar(r2, mut_res[:,1], color=my_cols[1], width=barWidth, edgecolor='white', label='A->C,T->G')
plt.bar(r3, mut_res[:,2], color=my_cols[2], width=barWidth, edgecolor='white', label='A->G,T->C')
plt.bar(r4, mut_res[:,3], color=my_cols[3], width=barWidth, edgecolor='white', label='G->T,C->A')
plt.bar(r5, mut_res[:,4], color=my_cols[4], width=barWidth, edgecolor='white', label='G->C,C->G')
plt.bar(r6, mut_res[:,5], color=my_cols[5], width=barWidth, edgecolor='white', label='C->T,G->A')
plt.ylabel('Frequency', fontweight='bold',fontsize=myfonts)
plt.xticks((np.array(r3)+np.array(r4))/2, ['Patient 04','Patient 09','Patient 12','Patient 15','Patient 16','Patient 26','Other Patients'],fontsize=myfonts)
plt.yticks(fontsize=myfonts)
plt.legend(fontsize=13)
plt.savefig('allS_w10SNP_allMutTypes_groupedByPatient.png',dpi=300,transparent=True)
plt.show()




##########################################################################################
##########################################################################################
######### Antimicrobial Resisistancies
######### Extended Data Figure 
##########################################################################################
##########################################################################################
# data generated /Users/u_key/Documents/mit/stapAD/readme/staphAD_AMR.sh l.1+

from matplotlib.backends.backend_pdf import PdfPages

########################################################
####### w Major Minor Lineage indicated
####### Basic QC reads for All Major/Minor Strains
####### use basic QC fastq
########################################################
# all isolates
d=pd.read_csv('/Users/u_key/Documents/mit/stapAD/mykrobe/allres_staphAD_basicQC.csv', index_col=0)

# minor lineage isolates 
minor_lineage_isolates=pd.read_csv('minor_lineage_isolates.txt',header=None)
minor_lineage_isolates[['subject', 'visit','location','isolate']] = minor_lineage_isolates[0].str.split('_', -1, expand=True)


fields = 5*['S','R','r','NA','S','R','r','NA','Major','Minor','None']
colors = 5*['#6EAF46','#045a8d', '#a6bddb','white','#6EAF46','#045a8d', '#a6bddb','white','black','lightgrey','white']
labels = 5*['S', 'R', 'r','NA','S', 'R', 'r','NA','Major','Minor','None']

subjects_numeric_sorted = sorted([int(i[1:]) for i in d['subject'].unique()])

with PdfPages('staphAD_basicQC_allMajorMinorIso_amr_per_visit_withMajMinor.pdf') as pdf:
    fig, ax = plt.subplots(0, figsize=(54, 32))
    n=0
    for px in subjects_numeric_sorted:
        p = 'S'+str(px)
        d_sub = d[d['subject']==p]
        perc_S_R_r_per_visit = np.array([])
        for ab in d_sub['drug'].unique():
            perc_maj_min_per_visit  = np.array([])
            for v in ['V1','V2','V3','V4','V5']:
                # record at each visit 3 possible states of susceptibility or NA (if visit had no data)
                if v in d_sub['visit'].unique():
                    res = np.array([])
                    res_major_minor_freq = np.array([])
                    d_sub_v_major_minor = d_sub[d_sub['visit']==v]
                    # get subset of strains 
                    d_sub_v = d_sub_v_major_minor[~d_sub_v_major_minor['isolate'].isin(minor_lineage_isolates['isolate'])] # major lineage isolates
                    d_sub_v_minor = d_sub_v_major_minor[d_sub_v_major_minor['isolate'].isin(minor_lineage_isolates['isolate'])] # minor lineage isolates
                    major_lineage_freq_visit = np.unique(d_sub_v['isolate']).size/np.unique(d_sub_v_major_minor['isolate']).size
                    minor_lineage_freq_visit = 1 - major_lineage_freq_visit
                    # get major lineage frequencies for stacked barplot
                    sus = np.unique(d_sub_v[d_sub_v['drug']==ab]['susceptibility'],return_counts=True)
                    for i,t in enumerate(['S','R','r','placeholder_for_NA_put_0']):
                        if t in sus[0]:
                            res = np.append(res,sus[1][np.where(sus[0]==t)[0]])
                        else:
                            res = np.append(res,0)
                    # get minor lineage frequencies for stacked barplot
                    sus = np.unique(d_sub_v_minor[d_sub_v_minor['drug']==ab]['susceptibility'],return_counts=True)
                    for i,t in enumerate(['S','R','r','placeholder_for_NA_put_0']):
                        if t in sus[0]:
                            res = np.append(res,sus[1][np.where(sus[0]==t)[0]])
                        else:
                            res = np.append(res,0)
                    res = np.append(res,[0,0,0]) # add 0 placeholder for minor/major/none freq
                    res_major_minor_freq = np.array([0,0,0,0,0,0,0,0,major_lineage_freq_visit,minor_lineage_freq_visit,0]) # store minor/major in separate vector
                else: # no data at visit -> NA = 1(00%)
                    res = np.array([0,0,0,1,0,0,0,0,0,0,0])  
                    res_major_minor_freq = np.array([0,0,0,0,0,0,0,0,0,0,1])
                res = res/np.sum(res) # get relative values
                perc_S_R_r_per_visit = np.concatenate((perc_S_R_r_per_visit,res))
                res_major_minor_freq = res_major_minor_freq/np.sum(res_major_minor_freq)
                perc_maj_min_per_visit = np.concatenate((perc_maj_min_per_visit,res_major_minor_freq))
        perc_S_R_r_per_visit = perc_S_R_r_per_visit.reshape(d_sub_v_major_minor['drug'].unique().size,55)
        perc_S_R_r_per_visit = np.vstack((perc_S_R_r_per_visit,perc_maj_min_per_visit)) # add major/minor to array with stacked info
        n += 1
        ax = fig.add_subplot(5,5,n)
        # fig, ax = plt.subplots(1, figsize=(12, 10))# plot bars
        plt.title('Patient '+str(px),loc='left',fontsize=20)
        left = (len(np.unique(d['drug']))+1) * [0]
        for idx, name in enumerate(fields):
            plt.barh(np.arange(len(np.unique(d['drug']))+1), perc_S_R_r_per_visit[:,idx], left = left, color=colors[idx])
            left = left + perc_S_R_r_per_visit[:,idx]
        # plt.legend(labels[0:4],bbox_to_anchor=([0.55, 1, 0, 0]), ncol=4, frameon=False)
        # if p == d['subject'].unique()[-1]:
        #     plt.legend(labels[0:4],fontsize=20)
        [plt.axvline(x=v, color="grey",linestyle="--") for v in [1,2,3,4]]
        ax.set_yticks(np.arange(13))
        ax.set_yticklabels(np.append(np.unique(d['drug']),'LINEAGE'),fontsize=14)
        plt.xlabel('Visit',fontsize=14)
        # Hide major tick labels
        ax.set_xticklabels('')
        # Customize minor tick labels
        ax.set_xticks([0.5,1.5,2.5,3.5,4.5], minor=True)
        ax.set_xticklabels(['1','2','3','4','5'], minor=True,fontsize=14)
    pdf.savefig(fig)




