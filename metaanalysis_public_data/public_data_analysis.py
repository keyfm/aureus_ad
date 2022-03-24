########################################################################################################################################################################################################
######################### Analysis of all 4 public datasets combined 
######################### 276 isolates 
########################################################################################################################################################################################################

####*** Prerequisite ***####
# - assemble and annotate each isolate (Extended Data Table 9 Key et al.) re-using the assembly snakemake presented for the within_host_analysis.
#   - samples should be appended by '_jac18'/'_bjd18'/'_jid18'/'elf17' to indicate the original study (tags used during downstream analysis)
# - use CL Blast (-outfmt 5 -max_hsps 1) and query all genes of interest (Extended Data Table 8) to each assembly
# --> The resulting genome.fasta, genome.gff, and query.xml files form the input for the code below


''' load libraries '''
import pandas as pd
import numpy as np
import argparse,sys,re
from Bio.Blast import NCBIXML
from Bio import SeqIO
import warnings, glob
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
import statsmodels.api as sm
import seaborn as sns
import scipy.stats as stats

''' 
    FUNCTIONS
    from functions_for_analysis.py
'''

SCRIPTS_DIRECTORY = "../modules"
sys.path.insert(0, SCRIPTS_DIRECTORY)

import analysispy_module as apy

# import importlib
# importlib.reload(apy)


''' 
    FUNCTIONS
    NEW
'''    

def subject_xml_paths(path):
    # search xml files: [path]/[subject]*[.xml]
    # return list for subject indentifiers. Assumed pattern: [path]/[subject_identifier].xml
    if path[-1] == "/":
        path = path[:-1] # remove trailing "/" if provided        
    file_paths_ls = glob.glob(path+'/*.xml')
    subject_xml_ls = [p.replace(path+"/", '') for p in file_paths_ls] # remove path string 
    subject_ls = [p.replace(".xml", '') for p in subject_xml_ls] # remove .xml strin >> subject-identifier list
    return [subject_ls,file_paths_ls]


def read_blast_records(file_paths_ls,xml_chr_style='prokka'):
    # loop over each subjects xml
    # report blast hit with highest bit score!
    blast_rec_per_subject_lopd = [] # list of pd.datframe records per subject. ordered as subject_ls!
    for subj_xml_path in file_paths_ls:
        print(subj_xml_path)
        with open(subj_xml_path) as xml:
            blast_records = NCBIXML.parse(xml)
            # loop over blast recored and add relevant info to dict (turned to pd.DF later)
            lod_blast_res = [] # add all blast_res dict for each record.
            for blast_record in blast_records:
                blast_res_dict = {} # store info per record
                blast_res_dict['query_fasta_header'] = blast_record.query.replace(",", "_") # replace possible commas (avoid interference with output csv)
                blast_res_dict['gene'] = blast_record.query.split('_')[0]
                blast_res_dict['ensemblID'] = '_'.join(blast_record.query.split(' ')[0].split('_')[1:]) # accounts for '_' w/in ensembl id
                blast_res_dict['length_ref'] = blast_record.query_letters # length of reference seq
                blast_res_dict['num_hits'] = len(blast_record.alignments)                
                # record blast hit(s) for query-seq
                if blast_res_dict['num_hits'] > 0:
                    lod_query_hits = []
                    for alignment in blast_record.alignments:
                        blast_query_hit = {}
                        if xml_chr_style == 'prokka':
                            blast_query_hit['chr_string'] = alignment.hit_id
                            blast_query_hit['chr_cov'] = alignment.hit_id.split('_')[-1] # add that to spot results from spurious assemblies
                            blast_query_hit['chr'] = alignment.hit_id.split('_')[1]
                        elif xml_chr_style == 'nakamuraSTM20':
                            blast_query_hit['chr_string'] = alignment.hit_id.split('|')[1]
                            blast_query_hit['chr_cov'] = np.nan # not reported in scaffold header by nakamuraSTM20
                            blast_query_hit['chr'] = int(alignment.hit_id.split('|')[1].split('.')[0][10:]) # extract part from string that contains scaffold number
                        for hsp in alignment.hsps:
                            blast_query_hit['score_e-value'] = hsp.expect
                            blast_query_hit['score_bit'] = hsp.bits
                            blast_query_hit['score_blast'] = hsp.score
                            blast_query_hit['loc1'] = hsp.sbjct_start
                            blast_query_hit['loc2'] = hsp.sbjct_end
                            blast_query_hit['length_query'] = len(hsp.query)
                            blast_query_hit['mismatches'] = hsp.match.count(' ')
                            if hsp.sbjct_start < hsp.sbjct_end:
                                blast_query_hit['strand'] = 1
                                blast_query_hit['gene_start'] = hsp.sbjct_start
                                blast_query_hit['gene_end'] = hsp.sbjct_end
                            else:
                                blast_query_hit['strand'] = -1
                                blast_query_hit['gene_start'] = hsp.sbjct_end
                                blast_query_hit['gene_end'] = hsp.sbjct_start
                        lod_query_hits.append(blast_query_hit)
                    # extract alignment with max score and add to blast_res_dict
                    blast_scores = [x['score_blast'] for x in lod_query_hits]
                    # print(bit_scores)
                    idx_best_alignment = blast_scores.index(max(blast_scores))
                    blast_res_dict.update(lod_query_hits[idx_best_alignment])
                # append results of xml to lod (> dataframe)
                lod_blast_res.append(blast_res_dict)                
        blast_record_pd = pd.DataFrame(lod_blast_res)
        blast_rec_per_subject_lopd.append(blast_record_pd)
    return blast_rec_per_subject_lopd




def match_blast2annotation(blast_table,scafNames,contig_lengths,annotation_genes,var_bp,min_percent_ref_aligned):
    # loop over supplied blast table (extracted from xml) and search for match in annotation based on chr/start/end coordinates    
    numeric_gene_status_ar = np.zeros( blast_table.shape[0] ,dtype=int ) # ls with element per gene in blast table. 0/1/2/3/4/5 
    # 0 noBlast
    # 1 Blast/noProkka
    # 2 Blast/Prokka
    # 3 Pseudogene (Blast hit spread across two different annotations)
    # 4 Blast start has annotated gene but gene end missing. possible truncated gene (pseudogene)
    # 5 as 4, but end of blast equals end of chr, thus assembly quality prevent explicit inference
    ugi = pd.Series(np.nan,index=range(blast_table.shape[0]),dtype=object) # pd series for inferred unified_gene_identifier_blastHit (if inexistent: NaN)
    for index, row in blast_table.iterrows():
        #print(index)
        # skip empty blast entry
        if row.isnull()['chr_string']:
            continue
        # skip if blast-located chr is empty annotation (prokka)
        chr_idx = [i for i,item in enumerate(scafNames) if row['chr_string'] in item] #get chr index; np.char.find() throws an error...
        # we accept blast hit, now search annotation for presence of gene
        if annotation_genes[chr_idx[0]].empty:
            numeric_gene_status_ar[index] = 1 # blast yes; prokka no
            continue
        # test if length of ref fasta was recovered in large by blast
        percent_ref_aligned = row['length_query']/row['length_ref']
        if percent_ref_aligned < min_percent_ref_aligned:
            # loss of query compared to ref due to assembly limitations?
            # <>3 bcs if gene at chr end, it stops at last complete triplett, thus chr end might be +-2b away
            if contig_lengths[ int(row['chr'])-1 ] > (row['gene_end']-3) and contig_lengths[ int(row['chr'])-1 ] < (row['gene_end']+3):
                numeric_gene_status_ar[index] = 5 # extra category: gene end == chr end....assembly quality prevent explicit inference: TruncationAtChrEnd
                continue
            else:
                continue # blast alignment covers insufficient ref seq
        # test if blast entry is also annotated by prokka
        s = annotation_genes[chr_idx[0]]['loc1'].between(row['gene_start']-var_bp, row['gene_start']+var_bp, inclusive=True)
        idx_match_start = s[s].index # index in annotated_genes[chr] where blast hit start matches; has to be same as idx_match_end
        s = annotation_genes[chr_idx[0]]['loc2'].between(row['gene_end']-var_bp, row['gene_end']+var_bp, inclusive=True)
        idx_match_end = s[s].index # index in annotated_genes[chr] where blast hit end matches; has to be same as idx_match_start
        #print(row['gene'])
        #print(idx_match_start)
        #print(idx_match_end)        
        # gene annotated at blast start pos
        if idx_match_start.size > 0:
            # with large var_bp values it can happen that multiple genes are hit.
            # > test if one id present in both > true > consider hit (I call that greedy)
            num_overlapping_ids = len(idx_match_start[np.in1d(idx_match_start,idx_match_end)])
            if num_overlapping_ids == 1:
                # start/end of blast have start/end of gene > blast hit also annotated
                numeric_gene_status_ar[index] = 2
                ugi[index] = annotation_genes[chr_idx[0]]['orthologtag'][ idx_match_start[np.in1d(idx_match_start,idx_match_end)] ].values[0] # add ortholog tag to list.
            elif num_overlapping_ids == 0 and idx_match_end.size > 0:
                # found annotated gene at start/end, but those are different genes > pseudogene. premature stop codon, sometimes leads prokka to annotate gene tail again (if "start codon" available)
                numeric_gene_status_ar[index] = 3
            elif idx_match_end.size == 0:
                # premature stop truncation but no new annotation of tail. Likely PG, but probably worth some beta-testing > assign for now own color!
                # <>3 bcs if gene at chr end, it stops at last complete triplett, thus chr end might be +-2b away
                if contig_lengths[ int(row['chr'])-1 ] > (row['gene_end']-3) and contig_lengths[ int(row['chr'])-1 ] < (row['gene_end']+3):
                    numeric_gene_status_ar[index] = 5 # extra category: gene end == chr end....assembly quality prevent explicit inference
                else:
                    numeric_gene_status_ar[index] = 4 # gene truncation. possible PG
            else:
                # does that ever happen?
                warnings.warn("Blast ok, but none of my annotation conditions fullfiled. Consider to lower var_bp.")
        else:
            numeric_gene_status_ar[index] = 1 # blast yes; prokka no
        blast_table['UnifiedGeneIdent'] = ugi # update blast table
    return [numeric_gene_status_ar,blast_table]

                
def heatmap_blast_anno_match(goi_blast_anno_pd,out_pdf,x_axis_label=None,horizontal_line_at=None,vertical_line_at=None,plot_title='nonspecified'):
    # generate heatmap for blast, blast-prokka presence of genes
    rc={'figure.figsize':(50,60), 'font.size': 32, 'axes.labelsize': 32, 'legend.fontsize': 32.0, 
    'axes.titlesize': 32, 'xtick.labelsize': 32, 'ytick.labelsize': 15}
    plt.rcParams.update(**rc)
    mycol = ('gray','skyblue','blue','orange','blueviolet','olive')
    ax = sns.heatmap(goi_blast_anno_pd, cmap=ListedColormap(mycol), linewidths=0.5, linecolor='lightgray',xticklabels=x_axis_label)
    ax.axes.set_title(plot_title)
    # Manually specify colorbar labelling after it's been generated
    colorbar = ax.collections[0].colorbar
    llv = (5/6)/2
    llvf = (5/6)
    colorbar.set_ticks([llv,llvf+llv,2*llvf+llv,3*llvf+llv,4*llvf+llv,5*llvf+llv])
    colorbar.set_ticklabels(['noBlast-noProkka', 'Blast-noProkka', 'Blast-Prokka','Pseudogene','Truncation','TruncationAtChrEnd']) # size 
    ax.figure.axes[-1].yaxis.label.set_size(50)
    if not np.any(horizontal_line_at == None):
        ax.hlines([horizontal_line_at], *ax.get_xlim(),colors='white',linewidth=10.0) # vertical line
    if not np.any(vertical_line_at == None):
        for i in vertical_line_at:
            ax.axvline(i, linewidth=5, c='w')
    fig = ax.get_figure()
#    plt.ylabel('AD                           NC')
    fig.savefig(out_pdf) 


def pattern_search(ar,pattern):
    ''' Returns index where pattern match in array '''
    indices = [i for i,x in enumerate(ar) if re.search(pattern,x)]
    return np.array(indices)
    # return ar.take(indices)

  


''' MAIN '''

if __name__ == "__main__":    
    # Specify paths and files
    xml_path = "" # Path to folder with all Blast results (in xml format)
    orth_anno_file = "annotation_orthologs.tsv" # file summarizes all orthologs across each isolate assembly. generated during snakemake-assembly
    ref_genome_folder = "" # folder to each isolate assembled genome file
    accepted_diff_bp = 100
    out_pdf = "heatmap_target_gene_presence_absence.pdf"

    # fix paths [xml path doen inside function]
    if ref_genome_folder[-1] == "/":
        ref_genome_folder = ref_genome_folder[:-1]

    # read ortholog inference of prokka annotations [unified-gene-identifier ...ugi...used as row index]
        with open(orth_anno_file) as f:
        ortholog_df = pd.read_table(orth_anno_file,index_col=0)
    
    ## parse blast xml records for each isolate
    # get list of subject_identifier and list of complete xml filepaths
    [subject_ls,xml_path_ls] = subject_xml_paths(xml_path)
        
    # read blast data from xml files    
    blast_rec_per_subject_lopd = read_blast_records(xml_path_ls)
        
    
    ####################################
    ### Infer agr & cap type and clean-up tables
    ####################################
    # clean up blast records: infer agr type and cap type based on max bit_score and remove non-existent agr/cap type genes 
    # use max bit score sum of agr / capHIJK to infer type
    # ATTENTION:  analysis specific to gene order in fasta file with query seq (SOM Table 8 in Key et al)
    agr_type_rows_lol = [[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]] # type i,ii,iii,iv order in capsule_goi_ensembl_nt_cap_trfak_list3
    cap58hijk_rows_lol = [[23,24,25,26],[32,33,34,35]] # cap5hijk, cap8hijk in capsule_goi_ensembl_nt_cap_trfak_list3
    blast_rec_per_subject_bestAgrCap_lopd = []
    subject_cap_agr_type_lol_elf17 = []
    for blast_rec in blast_rec_per_subject_lopd:
        ## infer agr/cap type based on best bit score
        bitS_max_agr_cap = []
        bit_score_sum = []
        for row_ls in agr_type_rows_lol:
            bit_score_sum.append(np.sum(blast_rec['score_bit'][row_ls]))
        bitS_max_agr_cap.append(bit_score_sum.index(max(bit_score_sum))) # get index of agr_type_rows_lol where best alignemnt agr is 
        bit_score_sum = []
        for row_ls in cap58hijk_rows_lol:
            bit_score_sum.append(np.sum(blast_rec['score_bit'][row_ls]))
        bitS_max_agr_cap.append(bit_score_sum.index(max(bit_score_sum))) # get index of cap58hijk_rows_lol where best alignemnt agr is 
        subject_cap_agr_type_lol_elf17.append( [['i','ii','iii','iv'][bitS_max_agr_cap[0]],['5','8'][bitS_max_agr_cap[1]]] )
        ## filter blast table for best agr/cap
        rows_for_removal = []
        for i in range(len(agr_type_rows_lol)):
            if not i==bitS_max_agr_cap[0]: # [0] agr index
                rows_for_removal.append(agr_type_rows_lol[i])
        for i in range(len(cap58hijk_rows_lol)):
            if not i==bitS_max_agr_cap[1]: # [1] cap index
                rows_for_removal.append(cap58hijk_rows_lol[i])
        rows_for_removal = [k for j in rows_for_removal for k in j] # lol to list
        clean_blast_rec = blast_rec.drop(rows_for_removal)
        ## reorder, so cap8hijk are in-frame
        clean_blast_rec = clean_blast_rec.reset_index(drop=True)
        if bitS_max_agr_cap[1] == 1: # ie. is cap8
            clean_blast_rec = clean_blast_rec.reindex([0,1,2,3,4,5,6,7,8,9,10,16,17,18,19,11,12,13,14,15,20,21,22,23,24,25,26,27,28,29,30,31,32]) # required capsule_goi_ensembl_nt_cap_trfak_list3!!
            clean_blast_rec = clean_blast_rec.reset_index(drop=True)
        blast_rec_per_subject_bestAgrCap_lopd.append( clean_blast_rec )
        
    
    ####################################
    ### Infer ALL blast hits in assembled genome-annotation for each subject
    ####################################
    # Here we match the actual annotation (by prokka) with the blast result. For each gene of interest the result is encoded in a single numeric:
    # 0 noBlast hit
    # 1 Blast hit but no Prokka annotation in that location
    # 2 Blast hit and Prokka annotation >> gene is present!
    # 3 Pseudogene (Blast hit spread across two annotated genes)
    # 4 Blast start has annotated gene but gene end missing. possible truncated gene (pseudogene)
    # 5 as 4, but end of blast equals end of chr, thus assembly quality prevent explicit inference

    goi_blast_anno_pd_clean_oth3elf17 = pd.DataFrame()
    for s,sub in enumerate(subject_ls):
        print(sub)
        ## read in genome information [NOTE: this is recycled from analysis.py in order to flawlessly integrate]        
        [chrStarts, genomeLength, scafNames] = apy.genomestats(ref_genome_folder+"/"+sub);
        contig_lengths = np.append(chrStarts[1:], genomeLength) - chrStarts # required to test if gene end coincides with chr end (suggesting incomplete assembly affect gene annotation)
        annotation_genes = apy.parse_gff(ref_genome_folder+"/"+sub,scafNames,ortholog_df[sub],forceReDo=False) # ref_genome_folder+"/annotation_genes.pandas.py.pk1"
        # estimate assembly quality based on mean chr coverage (assume each chr has unique coverage)
        print(pd.to_numeric(blast_rec_per_subject_bestAgrCap_lopd[s]['chr_cov'].dropna().unique()).mean())
        if pd.to_numeric(blast_rec_per_subject_bestAgrCap_lopd[s]['chr_cov'].dropna().unique()).mean() > 4:
            # loop over blast hits and search for annotated gene based on start/end coordinates
            [goi_blast_anno_pd_clean_oth3elf17[sub],blast_rec_per_subject_bestAgrCap_lopd[s]] = match_blast2annotation(blast_rec_per_subject_bestAgrCap_lopd[s],scafNames,contig_lengths,annotation_genes,var_bp=accepted_diff_bp,min_percent_ref_aligned=0.9)

    # format
    goi_blast_anno_pd_clean_oth3elf17 = goi_blast_anno_pd_clean_oth3elf17.transpose()
    goi_blast_anno_pd_sort2plot = goi_blast_anno_pd_clean_oth3elf17.sort_values(by=[3,7], ascending=[False,False]) # goi_blast_anno_wStatus_pd['Case_status_'].value_counts() >> 49x '0' --> horizontal line position
    
    # plot
    out_pdf = '/Users/u_key/Documents/mit/stapAD/pdf/public_data/sglIsolatePerSubject_PDonly/young105spl_elife17_match_blast2gff_ensembl_nt_cap_trfak_list3_cleanAgrCap.pdf'
    heatmap_blast_anno_match(goi_blast_anno_pd_sort2plot,out_pdf,x_axis_label=blast_rec_per_subject_bestAgrCap_lopd[0]['gene'],vertical_line_at=[4,20],plot_title='Young elife17')

    ## for downstram analysis
    blast_rec_per_subject_bestAgrCap_oth3elf17_lopd = blast_rec_per_subject_bestAgrCap_lopd

    


    ############################################################################################################
    ################ Analyses GOI and iTol for capD tree (main) and agr tree (SOM)
    ############################################################################################################
    ####################################
    ### Merge and final data prep
    ####################################
    ## merge elife and other3 datasets in order
    goi_blast_anno_pd_clean_oth3elf17 = pd.concat([goi_blast_anno_pd_clean_oth3, goi_blast_anno_pd_clean_elf17])    
    blast_rec_per_subject_bestAgrCap_oth3elf17_lopd = blast_rec_per_subject_bestAgrCap_oth3_lopd + blast_rec_per_subject_bestAgrCap_elf17_lopd    
    subject_cap_agr_type_lol_oth3elf17 = subject_ls
 
    ## drop rsaA rna gene
    goi_blast_anno_pd_clean_oth3elf17 = goi_blast_anno_pd_clean_oth3elf17.drop(columns=27)
    gene_ids = [gene for i,gene in enumerate(list(blast_rec_per_subject_bestAgrCap_oth3elf17_lopd[0]['gene'])) if i !=27 ] # remove rsaA here too
    # remove agr type specifier and cap 5/8 from gene list as this is not part of the information any longer
    gene_ids = [s.replace('8', '') for s in gene_ids] 
    gene_ids = [s.replace('5', '') for s in gene_ids]
    gene_ids = [s.replace('-i', '') for s in gene_ids]
    
    ############
    ### Build indices for each category
    ############
    # AD: 'AD*bjd18' 'PSAE*jid18'
    # OtherInf: *elife
    # healthy: 'NC*bjd18' .*NC.-jid18
    disease_status_all = np.zeros(subject_ls_oth3elf17.shape, dtype=object) # record for som table
    
    ## AD: 'P.*jid18|AD_.*bjd18' > 110
    p = re.compile('P.*jid18|AD_.*bjd18|.*jac18')
    ad_indices = np.array([ s for i,s in enumerate(goi_blast_anno_pd_clean_oth3elf17.index) if  p.match(s) ])
    disease_status_all[ np.array([ i for i,s in enumerate(goi_blast_anno_pd_clean_oth3elf17.index) if  p.match(s) ],dtype=int) ] = 'Atopic Dermatitis'
    
    ## otherInf: '.*blood|.*softtissue|.*bone-joint' > 99
    p = re.compile('.*blood|.*softtissue|.*bone-joint')
    oi_indices = np.array([ s for i,s in enumerate(goi_blast_anno_pd_clean_oth3elf17.index) if  p.match(s) ])
    disease_status_all[ np.array([ i for i,s in enumerate(goi_blast_anno_pd_clean_oth3elf17.index) if  p.match(s) ],dtype=int) ] = 'Other Infections'
    
    ## healthCtrl: 'NC.*bjd18|.*NC.-jid18' >> 67
    p = re.compile('NC.*bjd18|.*NC.-jid18')
    ct_indices = np.array([ s for i,s in enumerate(goi_blast_anno_pd_clean_oth3elf17.index) if  p.match(s) ])
    disease_status_all[ np.array([ i for i,s in enumerate(goi_blast_anno_pd_clean_oth3elf17.index) if  p.match(s) ],dtype=int) ] = 'Healthy'
 

    ####################################
    ### SOM Table public data results
    ####################################
    # Table contains: study,sample ID, cap/agr type; status GOI (0:absent, 1:present)

    ## build series with pub source
    isolate_source = np.zeros(subject_ls_oth3elf17.shape,dtype=object)
    s = pd.Series(subject_ls_oth3elf17)
    isolate_source[s.str.endswith('jac18')] = 'Edslev2018'
    isolate_source[s.str.endswith('bjd18')] = 'Harkins2018a'
    isolate_source[s.str.endswith('jid18')] = 'Harkins2018b'
    isolate_source[s.str.endswith('elf17')] = 'Young2017'

    ## build series with sample id
    p = re.compile('-jid18|-bjd18|-jac18|-elf17')
    sample_ids = s.str.replace(p, '', regex=True)

    som_table_arr = np.array([],dtype=object)
    row_ctr = 0
    for index, row in goi_blast_anno_pd_clean_oth3elf17.iterrows():
        goi_blast_simple = np.array([1 if j==2 else 0 for j in row],dtype=object) # put blast results to 0(absent) and 1(present)       
        if row_ctr == 0:
            som_table_arr = np.concatenate([[isolate_source[row_ctr]],[sample_ids[row_ctr]],[disease_status_all[row_ctr]],np.array(subject_cap_agr_type_lol_oth3elf17[row_ctr]),goi_blast_simple])
        else:
            som_table_arr = np.vstack( (som_table_arr,np.concatenate([[isolate_source[row_ctr]],[sample_ids[row_ctr]],[disease_status_all[row_ctr]],np.array(subject_cap_agr_type_lol_oth3elf17[row_ctr]),goi_blast_simple]) ) )
        row_ctr += 1
    
    
    ## sort 
    som_table_df = pd.DataFrame(som_table_arr,columns=['Source','Isolate ID','Disease Status','agr type','capsule type']+gene_ids ).sort_values(['Source', 'Isolate ID'], ascending = [True, True])
    with open('/Users/u_key/Documents/mit/stapAD/tables/public_data/public_data_276_results.csv','w') as file:
        som_table_df.to_csv(file,header=True,index=False)


    ####################################
    ### PLOTTING
    ####################################
    # Heatmap: % non-func gene in each category (AD,healthy,otherind)
    # remove RNA gene rsaA
    # stacked bar horizontal: # of isolates per pd dataset (not used at the end in ms)
    # barplot with %-non.functional capD and agrC incl. stats and error bars
    ## assemble three plots in illustrator
   

    #########
    ### heatmap all GOI
    #########
    ## build np matrix with % non func per category
    # == 2 is functional gene inferred. All other are not!    
    perc_nonFunc_gene = np.sum(np.array(goi_blast_anno_pd_clean_oth3elf17.loc[oi_indices] != 2),axis=0)/len(oi_indices)
    perc_nonFunc_gene = np.vstack((perc_nonFunc_gene, np.sum(np.array(goi_blast_anno_pd_clean_oth3elf17.loc[ct_indices] != 2),axis=0)/len(ct_indices)))
    perc_nonFunc_gene = np.vstack((perc_nonFunc_gene,np.sum(np.array(goi_blast_anno_pd_clean_oth3elf17.loc[ad_indices] != 2),axis=0)/len(ad_indices)))

    # get order only for heatmap: cap,agr,other
    new_order = [ 4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 0,1,2,3,20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31] # cap,agr, other genes
    np.array(gene_ids)[new_order] # sanity check
    
    ## plotting
    fig, ax = plt.subplots(figsize=(10,3))
    plt.rc('font', family='Helvetica');plt.rcParams['pdf.fonttype'] = 42 # change font to helvetica
    im = ax.imshow(perc_nonFunc_gene[:,np.array(new_order)],cmap='Blues')

    # We want to show all ticks...
    ax.set_xticks(np.arange(perc_nonFunc_gene.shape[1]))
    ax.set_yticks(np.arange(perc_nonFunc_gene.shape[0]))
    # ... and label them with the respective list entries
    ax.set_xticklabels(np.array(gene_ids)[new_order])
    ax.set_yticklabels(['OI','CT','AD'])
    plt.xticks(fontsize=14 )
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90) #, ha="right",
             # rotation_mode="anchor")
    plt.colorbar(im)
    plt.savefig('components_multipanel_fig/heatmap_v1.svg',transparent=True)
    
    ### plot # of isolates per pd dataset in brown/green/hollow coloring as in tree
    # NOT USED in PUB
    count_ad_oi_ct = np.array([len(ad_indices),len(oi_indices),len(ct_indices)])

    plt.style.use('default')
    plt.axis('off')
    fig, ax = plt.subplots(figsize=(4, 10))
    # plt.bar(pos, prop_capd_pos, width=0.8, label='capD-', color='darkgray', bottom=prop_capd_neg)
    barlist=plt.bar(np.arange(3), count_ad_oi_ct, width=0.8, color=['#018571','#a6611a','white'], edgecolor=['#018571','#a6611a','#018571'],linewidth=3)
    # plt.ylim(0, 1)
    plt.xticks(np.arange(3), np.array(['Atopic\nDermatitis','Other\nInfection','Healthy']))
    plt.ylabel("Number isolates",fontsize=13)
    ax.tick_params(axis='both', which='major', labelsize=13)
    fig.subplots_adjust(bottom=0.25) 
    plt.setp(plt.gca().get_xticklabels(), rotation=0, horizontalalignment='center') # rotate axis labels
    plt.savefig('components_multipanel_fig/barplot_isolateCount_v1.pdf',transparent=True)

    ### plot gene-loss count all GOI (bargraph)
    # NOT USED in PUB
    # subset by category, color as in F3
    ad_count = np.sum(np.array(goi_blast_anno_pd_clean_oth3elf17.loc[ad_indices] != 2),axis=0)
    ct_count = np.sum(np.array(goi_blast_anno_pd_clean_oth3elf17.loc[ct_indices] != 2),axis=0)
    oi_count = np.sum(np.array(goi_blast_anno_pd_clean_oth3elf17.loc[oi_indices] != 2),axis=0)
    
    # matrices for later
    ad_ct_oi_count_neg = np.vstack((ad_count,ct_count,oi_count))
    ad_ct_oi_count_pos = np.vstack((len(ad_indices)-ad_count,len(ct_indices)-ct_count,len(oi_indices)-oi_count))

    ## plot
    width = 0.8
    fig, ax = plt.subplots(figsize=(10,3))
    x_bars = np.arange(goi_blast_anno_pd_clean_oth3elf17.shape[1])
    p1 = plt.bar(x_bars,oi_count,width,color='#a6611a', edgecolor='#a6611a')
    p2 = plt.bar(x_bars,ct_count,width,bottom=oi_count,color='white', edgecolor='#018571')
    p3 = plt.bar(x_bars,ad_count,width,bottom=(ct_count+oi_count),color='#018571', edgecolor='#018571')
    
    plt.ylabel('Number isolates with gene absent')
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
    plt.legend((p3[0], p2[0], p1[0]), ('Atopic Dermatitis', 'Healthy','Other Infection'))
    plt.savefig('barplot_noGeneCount_v2.svg',transparent=True)
    

    ################################################################################################
    ###### BARPLOT Fig 4
    ###### - %capd_negative in AD/controls/other infections
    ###### - %agrc_negative in AD/controls/other infections
    ################################################################################################
    # numbers:
    import statsmodels.api

    #### capD    
    ## get data
    # bar shows capD- capD+ proportion 
    # for AD, CTRL, other (elife)
    capd_idx = np.where(np.array(gene_ids)=='capD')[0] # [0] to remove tuple struc
    capd_neg_freq = np.flip(perc_nonFunc_gene[:,capd_idx]) # order AD/CT/OI
    capd_neg_count = ad_ct_oi_count_neg[:,capd_idx] # order AD/CT/OI
    capd_pos_count = ad_ct_oi_count_pos[:,capd_idx] # order AD/CT/OI
    
    pos = np.arange(3)
    
    # get sampling proportion for CI error bars
    ci = statsmodels.stats.proportion.proportion_confint(count=capd_neg_count[:,0],nobs=(capd_neg_count[:,0]+capd_pos_count[:,0]),alpha=0.05,method='normal')
    # build 2d array with error bar + - from prop_capd_neg value
    ci=np.vstack(ci)
    ci[0]=capd_neg_freq[:,0]-ci[0]
    ci[1]=ci[1]-capd_neg_freq[:,0]
    

    fig, ax = plt.subplots(figsize=(6, 4))
    plt.rc('font', family='Helvetica');plt.rcParams['pdf.fonttype'] = 42 # change font to helvetica
    # plt.bar(pos, prop_capd_pos, width=0.8, label='capD-', color='darkgray', bottom=prop_capd_neg)
    barlist=plt.bar(pos, capd_neg_freq[:,0], width=0.8, label='capD+', color=['#018571','white','#a6611a'], yerr=ci,edgecolor=['#018571','#018571','#a6611a'],linewidth=3)
    plt.hlines(0.5,0,1,color='black')
    plt.hlines(0.7,0,2,color='black')
    # plt.axhline(y=0.24,color='grey',linestyle='dashed') # capD mut in staphAD: 0.24: 6/25
    plt.ylim(0, 1)
    plt.xticks(pos, np.array(['Atopic\nDermatitis','Healthy','Other\nInfection']))
    plt.ylabel("Percent capD truncation",fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.subplots_adjust(bottom=0.25) 
    plt.setp(plt.gca().get_xticklabels(), rotation=0, horizontalalignment='center') # rotate axis labels
    plt.savefig('public_data_4datasets_capD_status.pdf',transparent=True)
    plt.show()

    #### agrC    
    ## get data
    # bar shows agrc- agrc+ proportion 
    # for AD, CTRL, other (elife)
    agrc_idx = np.where(np.array(gene_ids)=='agrC')[0] # [0] to remove tuple struc
    agrc_neg_freq = np.flip(perc_nonFunc_gene[:,agrc_idx]) # order AD/CT/OI
    agrc_neg_count = ad_ct_oi_count_neg[:,agrc_idx] # order AD/CT/OI
    agrc_pos_count = ad_ct_oi_count_pos[:,agrc_idx] # order AD/CT/OI
    
    # get sampling proportion for CI error bars
    ci = statsmodels.stats.proportion.proportion_confint(count=agrc_neg_count[:,0],nobs=(agrc_neg_count[:,0]+agrc_pos_count[:,0]),alpha=0.05,method='normal')
    # build 2d array with error bar + - from prop_agrc_neg value
    ci=np.vstack(ci)
    ci[0]=agrc_neg_freq[:,0]-ci[0]
    ci[1]=ci[1]-agrc_neg_freq[:,0]
    

    fig, ax = plt.subplots(figsize=(6, 4))
    plt.rc('font', family='Helvetica');plt.rcParams['pdf.fonttype'] = 42 # change font to helvetica
    # plt.bar(pos, prop_agrc_pos, width=0.8, label='agrc-', color='darkgray', bottom=prop_agrc_neg)
    barlist=plt.bar(np.arange(len(agrc_neg_freq)), agrc_neg_freq[:,0], width=0.8, label='agrc+', color=['#018571','white','#a6611a'], yerr=ci,edgecolor=['#018571','#018571','#a6611a'],linewidth=3)
    plt.hlines(0.35,0,1,color='black')
    plt.hlines(0.55,0,2,color='black')
    # plt.axhline(y=0.24,color='grey',linestyle='dashed') # agrc mut in staphAD: 0.24: 6/25
    plt.ylim(0, 1)
    plt.xticks(pos, np.array(['Atopic\nDermatitis','Healthy','Other\nInfection']))
    plt.ylabel("Percent agrC truncation",fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.subplots_adjust(bottom=0.25) 
    plt.setp(plt.gca().get_xticklabels(), rotation=0, horizontalalignment='center') # rotate axis labels
    plt.savefig('public_data_4datasets_agrc_status.pdf',transparent=True)
    plt.show()

    ##########################################################################################
    ##########################################################################################
    ######### Barplot public data overview
    ######### Figure 4a
    ##########################################################################################
    ##########################################################################################
    ### get data
    ## 48 AD Harkin 2018a BJD
    p = re.compile('AD_.*bjd18')
    np.array([ s for i,s in enumerate(goi_blast_anno_pd_clean_oth3elf17.index) if  p.match(s) ]).shape 

    ## 49 CT Harkin 2018a BJD
    p = re.compile('NC.*bjd18')
    np.array([ s for i,s in enumerate(goi_blast_anno_pd_clean_oth3elf17.index) if  p.match(s) ]).shape 

    ## 9 AD Harkin 2018b JID
    p = re.compile('P.*jid18')
    np.array([ s for i,s in enumerate(goi_blast_anno_pd_clean_oth3elf17.index) if  p.match(s) ]).shape 

    ## 18 CT Harkin 2018b JID
    p = re.compile('.*NC.-jid18')
    np.array([ s for i,s in enumerate(goi_blast_anno_pd_clean_oth3elf17.index) if  p.match(s) ]).shape 

    ## 53 AD Edslev 2018 JAC
    p = re.compile('.*jac18')
    np.array([ s for i,s in enumerate(goi_blast_anno_pd_clean_oth3elf17.index) if  p.match(s) ]).shape 

    ## 99 OI Young 2017 elife
    p = re.compile('.*blood|.*softtissue|.*bone-joint')
    np.array([ s for i,s in enumerate(goi_blast_anno_pd_clean_oth3elf17.index) if  p.match(s) ]).shape 


    # data arrays    
    sick = np.array([99,53,9,48]) #elife17, 
    heal = np.array([18,49])

    
    barwidth = 0.35
    xsick = [x+barwidth for x in np.arange(4)]
    xheal = [x for x in np.arange(2,4)]
    fig, ax = plt.subplots(figsize=(6, 4))
    plt.rc('font', family='Helvetica');plt.rcParams['pdf.fonttype'] = 42 # change font to helvetica
    plt.rcParams['axes.linewidth'] = 2 #set the value globally
    plt.barh( [xsick[x] for x in np.arange(2,4)], [sick[x] for x in np.arange(2,4)], color=['#018571','#018571','#018571'], height=barwidth, edgecolor='white', label='Atopic Dermatitis',)
    plt.barh(xheal, heal, color='white', height=barwidth, edgecolor="#018571", label='Healthy',linewidth=2)
    plt.barh(xsick[1]-barwidth/2, sick[1], color=['#018571'], height=barwidth, edgecolor='white')
    plt.barh(xsick[0]-barwidth, sick[0], color=['#a6611a'], height=barwidth, edgecolor='white', label='Other Infection',)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Add yticks on the middle of the group bars
    ylabelsTick = [x + (barwidth/2) for x in range(len(xsick))]
    plt.ylim(-0.5,4)
    ylabelsTick[0] = ylabelsTick[0] - (barwidth/2) # Young
    plt.yticks(ylabelsTick, ['Young 2017', 'Edslev 2018', 'Harkins 2018b', 'Harkin 2018a'],fontsize=16)
    plt.xticks(fontsize=16)
    plt.legend(prop={'size': 14})
    # plt.rcParams.update({'font.size': 122})
    plt.savefig('barplot_public data.pdf',bbox_inches='tight')



    
    ####################
    #### Fishers test on agrC and CapD
    ####################
    ### capD
    # ad loss capD: 25 from 110 isolates
    # ct loss capD: 1 from 67 isolates
    # oi loss capD: 11 from 99 isolates

    ## fisher: AD vs Healthy
    oddsratio, pvalue = stats.fisher_exact([[25, 1], [110-25, 67-1]],alternative='greater') 
    # (19.41176470588235, 2.536867391419118e-05)

    ## fisher AD vs. other Inf
    oddsratio, pvalue = stats.fisher_exact([[25, 11], [110-25, 99-11]],alternative='greater') 
    # (2.3529411764705883, 0.019947753123933258)

    ### agrC
    # ad loss agrC: 11 from 110 isolates
    # ct loss agrC: 2 from 67 isolates
    # oi loss agrC: 0 from 99 isolates

    import scipy.stats as stats
    ## fisher: AD vs Healthy
    oddsratio, pvalue = stats.fisher_exact([[11, 2], [110-11, 67-2]],alternative='greater')
    # (3.611111111111111, 0.07040517416928116) # two-sided: (3.611111111111111, 0.13509660497937007)

    ## fisher AD vs. other Inf
    oddsratio, pvalue = stats.fisher_exact([[11, 0], [110-11, 99-0]],alternative='greater')
    # (inf, 0.0006692004709397211) # two-sided: (inf, 0.0008670534658342698)

    ### agr-operon
    # ad loss agrC: 14 from 110 isolates
    # ct loss agrC: 3 from 67 isolates
    # oi loss agrC: 2 from 99 isolates

    import scipy.stats as stats
    ## fisher: AD vs Healthy
    oddsratio, pvalue = stats.fisher_exact([[14, 3], [110-14, 67-3]],alternative='greater')
    # (3.111111111111111, 0.057001948736045444) # two sided: (3.111111111111111, 0.11228589821540903)

    ## fisher AD vs. other Inf
    oddsratio, pvalue = stats.fisher_exact([[14, 2], [110-14, 99-2]],alternative='greater')
    # (7.072916666666667, 0.0028861589385491943) # two-sided: (7.072916666666667, 0.003604824562824347)


    ## test all genes for enrichment AD vs. CT or OI
    # P<0.001
    for i,ad_num in enumerate(ad_count):
        oddsratio, pvalue = stats.fisher_exact([[ad_num, ct_count[i]], [len(ad_indices)-ad_num, len(ct_indices)-ct_count[i]]],alternative='greater') 
        if pvalue < 0.001:
            print(gene_ids[i],'AD vs. CT',pvalue)
        oddsratio, pvalue = stats.fisher_exact([[ad_num, oi_count[i]], [len(ad_indices)-ad_num, len(oi_indices)-oi_count[i]]],alternative='greater')
        if pvalue < 0.001:
            print(gene_ids[i],'AD vs. OI',pvalue)

    # RESULTS
    # agrC AD vs. OI 0.0006692004709397211
    # capD AD vs. CT 0.00015752390873181136


    ############
    ### Fisher's test w/o Edslev-clonal outbreak study
    ############
    # AD: 'AD*bjd18' 'PSAE*jid18'
    # OtherInf: *elife
    # healthy: 'NC*bjd18' .*NC.-jid18
    
    ## Remove Edslev from data, because isolates apparently represent a somewhat clonal outbreak of CC1 isolates in Danish AD patients
    # Edslev provided only AD data, thus only AD related numbers change. OI/CT remain unchanged.
    # AD: 'P.*jid18|AD_.*bjd18' > 111
    p = re.compile('P.*jid18|AD_.*bjd18')
    ad_indices_noEdslev = np.array([ s for i,s in enumerate(goi_blast_anno_pd_clean_oth3elf17.index) if  p.match(s) ])
    ad_count_noEdslev = np.sum(np.array(goi_blast_anno_pd_clean_oth3elf17.loc[ad_indices_noEdslev] != 2),axis=0)

    ## test all genes for enrichment AD vs. CT or OI
    for i,ad_num in enumerate(ad_count_noEdslev):
        oddsratio, pvalue = stats.fisher_exact([[ad_num, ct_count[i]], [len(ad_indices_noEdslev)-ad_num, len(ct_indices)-ct_count[i]]],alternative='greater')
        if pvalue < 0.05:
            print(gene_ids[i],'AD vs. CT',pvalue)
        oddsratio, pvalue = stats.fisher_exact([[ad_num, oi_count[i]], [len(ad_indices_noEdslev)-ad_num, len(oi_indices)-oi_count[i]]],alternative='greater')
        if pvalue < 0.05:
            print(gene_ids[i],'AD vs. OI',pvalue)
    ## Results:
    # agrC AD vs. OI 0.005802072690118247
    # capD AD vs. CT 0.013582137033883978

    ############
    ### Fisher's test w/o CC8 (USA300) clone
    ############
    # Reviewer argues clonal outbreak should be removed.
    # a. Remove CC8 representative of USA300 (however, not all cc1 are USA300!)
    # b. Remove CC8-capDnegative-only as more refined representation of USA300-only
    
    ### a. Remove CC8 representative of USA300 (however, not all cc8 are USA300!)
    ### capD (substract CC8 representative of USA300, Note: CC8 is not ONLY USA300)
    # ad loss capD: 16 (25-9) from 99 (110-11) isolates
    # ct loss capD: 0 (1-1) from 66 (67-1) isolates
    # oi loss capD: 5 (11-6) from 91 (99-8) isolates

    ## fisher: AD vs Healthy
    oddsratio, pvalue = stats.fisher_exact([[16, 0], [99-16, 66]],alternative='greater')
    # (inf, 0.00016610725161645643)

    ## fisher AD vs. other Inf
    oddsratio, pvalue = stats.fisher_exact([[16, 5], [99-16, 91-5]],alternative='greater')
    # (3.31566265060241, 0.015965889929804696)
    
    ### b. Remove CC8-capDnegative-only as more refined representation of USA300-only
    ### capD 
    # ad loss capD: 16 (25-9) from 101 (110-9) isolates
    # ct loss capD: 0 (1-1) from 66 (67-1) isolates
    # oi loss capD: 5 (11-6) from 93 (99-6) isolates

    ## fisher: AD vs Healthy
    oddsratio, pvalue = stats.fisher_exact([[16, 0], [101-16, 66]]) # default two-sided!
    # (inf, 0.0002422675320801615)

    ## fisher AD vs. other Inf
    oddsratio, pvalue = stats.fisher_exact([[16, 5], [101-16, 93-5]]) # default two-sided!
    # (3.3129411764705883, 0.021382116201545114)
    


    ################################################    
    ### Fishers any capsule gene loss (but capD) AD vs healthy
    ################################################
    
    ## entire capsule locus
    gene_ids[4:20] > capsule genes
    cap_del_ad = np.sum(np.sum(np.array(goi_blast_anno_pd_clean_oth3elf17.iloc[:,4:20].loc[ad_indices] != 2),axis=1) > 0) # 31
    cap_del_oi = np.sum(np.sum(np.array(goi_blast_anno_pd_clean_oth3elf17.iloc[:,4:20].loc[oi_indices] != 2),axis=1) > 0) # 13
    cap_del_ct = np.sum(np.sum(np.array(goi_blast_anno_pd_clean_oth3elf17.iloc[:,4:20].loc[ct_indices] != 2),axis=1) > 0) # 3
    
    ## fisher: AD vs Healthy
    oddsratio, pvalue = stats.fisher_exact([[31, 3], [110-31, 67-3]],alternative='greater')
    # (8.371308016877638, 3.771124269323601e-05)

    ## fisher: AD vs OI
    oddsratio, pvalue = stats.fisher_exact([[31, 13], [110-31, 99-13]],alternative='greater')
    # (2.5959104186952286, 0.005851002345969357)


    ## capsule locus BUT capD
    cap_butD_del_ad = np.sum(np.sum(np.array(goi_blast_anno_pd_clean_oth3elf17.drop(7,axis=1).iloc[:,4:19].loc[ad_indices] != 2),axis=1) > 0) # 19
    cap_butD_del_oi = np.sum(np.sum(np.array(goi_blast_anno_pd_clean_oth3elf17.drop(7,axis=1).iloc[:,4:19].loc[oi_indices] != 2),axis=1) > 0) # 4
    cap_butD_del_ct = np.sum(np.sum(np.array(goi_blast_anno_pd_clean_oth3elf17.drop(7,axis=1).iloc[:,4:19].loc[ct_indices] != 2),axis=1) > 0) # 2

    ## fisher: AD vs Healthy
    oddsratio, pvalue = stats.fisher_exact([[19, 2], [110-19, 67-2]],alternative='greater')
    # (6.785714285714286, 0.0026679899061262976)

    ## fisher: AD vs OI
    oddsratio, pvalue = stats.fisher_exact([[19, 4], [110-19, 99-4]],alternative='greater')
    # (4.958791208791209, 0.0017401554310356544)

    ################################################    
    ### BUILD data structure for itol
    ################################################
    # iTol: agr type color
    # iTol: agr loss: differentiate ABCD
    # iTol: public data source (same as capD tree)
    # iTol: category AD/CT/OI (same as capD tree)
    # goi_blast_anno_pd_all = pd.concat([goi_blast_anno_pd_elf17,goi_blast_anno_pd_ot3]) # OLD
    # goi_blast_anno_pd_clean_oth3elf17 = pd.concat([goi_blast_anno_pd_clean_oth3, goi_blast_anno_pd_clean_elf17]) # NEW. Created above
    
    
    
    ### colorstrip for publication: capD and agr Tree
    # LEGEND_TITLE,Dataset_legend
    # LEGEND_SHAPES,2,2,2,2
    # LEGEND_COLORS,#FFFFFF,#C0C0C0,#808080,#000000
    # LEGEND_LABELS,Edslev_JAC18,Young_eLife17,Harkins_BJD18,Harkins_JID18

    #color_label = {'elf17':'#C0C0C0', 'jac18':'#FFFFFF','bjd18':'#808080','jid18':'#000000'}
    color_label = {'elf17':'#969696', 'jac18':'#d9d9d9','bjd18':'#808080','jid18':'#000000'} # grey scale, no white
    popup_label = {'elf17':'Young_eLife17', 'jac18':'Edslev_JAC18','bjd18':'Harkins_BJD18','jid18':'Harkins_JID18'}
    
    for s in subject_ls_oth3elf17:
        print(','.join([s,color_label[s.split('-')[-1]],popup_label[s.split('-')[-1]]]))
    
    # one sample missing? > USA300 stays uncolored
    
    
    ### binary data category - capD/agr Tree 
    ## AD and capD: 
    # all but otherInfection
    # checked  isolates for cryptic capD status!!!!
    for s in subject_ls_oth3elf17:
        s_noElf17 = s.replace('-elf17','') # pandas dataframe different format for elf17 data
        capd_status = goi_blast_anno_pd_clean_oth3elf17.loc[s_noElf17,7] # Attention: index is 0-based position of capD info
        if capd_status == 2:
            capd = -1
        else:
            capd = 1
        if s.startswith('AD'):
            ad_status = 1
        elif s.startswith('NC'):
            ad_status = 0
        elif s.split('-')[-1] == 'elf17':
            ad_status = -1
        elif s.split('_')[-1].startswith('ls'):
            ad_status = 1
        elif s.split('-')[-1] == 'jid18' and s.startswith('P'):
            ad_status = 1
        else:
            ad_status = 0
        print(','.join([s,str(ad_status),str(capd)]))
        

    ### binary data: otherInfection
    for s in subject_ls_oth3elf17:
        if not s.split('-')[-1] == 'elf17':
            print(','.join([s,str(1)]))
        else:
            print(','.join([s,str(-1)]))
                



    ### iTol input for CC inference colorstrip    
    # mlst inference bbased in srst2 (data provided on github)
    # color scheme. added grey for CC==NA. source: https://colorbrewer2.org/#type=qualitative&scheme=Pastel1&n=8
    cc_color_dc = {'CC1':'#fbb4ae','CC30':'#b3cde3','CC15':'#ccebc5','CC45':'#decbe4','CC5':'#fed9a6','CC8':'#ffffcc','CC97':'#e5d8bd','CC22':'#fddaec','NA':'#ffffff'}
    with open('srst2_pd_sglIso_res_wCC.txt','r') as file:
        for i,line in enumerate(file):
            line = line.strip().split(',')
            if i != 0 and line[0] in subject_ls_oth3elf17: # skip header & removed isolates
                    # print(line[0],'range',cc_color_dc[line[13]],line[13],'1') # colored tree
                    print(line[0],cc_color_dc[line[13]],line[13]) # colorstrip
    
    ### iTol input for AGR type inference colorstrip    
    agr_color_dc = {'agrA-i':'#f0f9e8','agrA-ii':'#7bccc4','agrA-iii':'#2b8cbe','agrA-iv':'#253494'}    
    for i,spl in enumerate(subject_ls_oth3elf17):
        agr_id = blast_rec_per_subject_bestAgrCap_oth3elf17_lopd[i]['gene'][0] # blast record stored gene name which contains inferred agr type
        agr_type = 'agr-'+agr_id.split('-')[1]
        print(spl,agr_color_dc[agr_id],agr_type)

    ### agr-genes absence binary
    ## triangle for agr loss. 4 colors for A/B/C/D loss differentiation
        ### binary data category - capD/agr Tree 
    ## AD and capD: 
    # all but otherInfection
    # checked  isolates for cryptic capD status!!!!
    # colors: A,B,C,D: ['#fdbe85','#fd8d3c','#e6550d','#a63603']
    for s in subject_ls_oth3elf17:
        s_noElf17 = s.replace('-elf17','') # pandas dataframe different format for elf17 data
        agrA_status = goi_blast_anno_pd_clean_oth3elf17.loc[s_noElf17,0]
        agrB_status = goi_blast_anno_pd_clean_oth3elf17.loc[s_noElf17,1]
        agrC_status = goi_blast_anno_pd_clean_oth3elf17.loc[s_noElf17,2]
        agrD_status = goi_blast_anno_pd_clean_oth3elf17.loc[s_noElf17,3]
        if agrD_status != 2: # exchange agrX_status to get res for each iTol file!
            print(s+',1')


