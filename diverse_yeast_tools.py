import os
import numpy as np
import scipy.stats as stats
import pandas as pd
import json
from collections import Counter
import pickle
import warnings
from Bio import SeqIO
from Bio import AlignIO   #, Align
from ete3 import Tree, EvolTree

base_dir = os.path.normpath('G:/My Drive/Crick_LMS/projects/diverse_yeasts/alphafold')
y1000plus_dir = os.path.normpath('C:/Users/heineib/Documents/GitHub/y1000plus_tools')
genomes_dir = os.path.normpath('G:/My Drive/Crick_LMS/external_data/genomes')

ctg_clade_dict = {'candida_albicans': 'CUG-Ser1',
                  'debaryomyces_hansenii': 'CUG-Ser1',  
                  'candida_tropicalis': 'CUG-Ser1',
                  'ascoidea_rubescens': 'CUG-Ser2',
                  'pachysolen_tannophilus': 'CUG-Ala'
                 }

ctg_clade_alt_res = {'CUG-Ala': 'A', 
                     'CUG-Ser1': 'S',
                     'CUG-Ser2': 'S'
                     }


model_spec_lookup = {'Scer': 'saccharomyces_cerevisiae',
                     'Calb': 'candida_albicans',
                     'Spom': 'schizosaccharomyces_pombe'
}

#Functions to parse sequence id from first line of fasta file in either shen or uniprot proteomes. 
def gene_id_shen(seq_record):
    gene_id = seq_record.description.split()[1].split('=')[1]
    return gene_id
    
def gene_id_uniprot(seq_record):
    gene_id = seq_record.description.split()[0]
    return gene_id

gene_id_function_dict = {
    "shen": gene_id_shen,
    "uniprot": gene_id_uniprot
}

def gene_id_retrieve(study, seq_record):
    gene_id = gene_id_function_dict[study](seq_record)
    return(gene_id)

def read_SGD_features(sgd_features_fname=genomes_dir + os.sep + os.path.normpath('saccharomyces_cerevisiae/SGD_features.tab')):
    
    #SC_features_fname = os.path.normpath(data_processing_dir + "/ortholog_files_regev/SGD_features.tab")
    SC_features_fname = os.path.normpath(sgd_features_fname)

    SC_features = pd.read_csv(SC_features_fname, sep = '\t', header=None)
    SC_orfs = SC_features.groupby(1).get_group('ORF')
    
    #Makes a dictionary to look up orfs by gene names.  This won't include all orfs - those without names had NaN in column 4 so 
    #are presumably left out. 
    SC_orfs_lookup = dict(zip(SC_orfs[4], SC_orfs[3]))
    SC_orfs_lookup[np.nan] = np.nan
    SC_genename_lookup = dict(zip(SC_orfs[3], SC_orfs[4]))
    SC_features_lookup = dict(zip(SC_orfs[3], SC_orfs[15]))
       
    return SC_orfs_lookup, SC_genename_lookup, SC_features_lookup

def SC_common_name_lookup(gene_list):
    #SC Common Name lookup
    #Input is a list of orfs, output is a list of common names, or if there is no common name, the orf name
    
    SC_orfs_lookup, SC_common_name_lookup, SC_features_lookup = read_SGD_features()
    SC_common_names = []
    
    for gene in gene_list: 
        try:  
            SC_common_name = SC_common_name_lookup[gene]
            if isinstance(SC_common_name,float):
                if math.isnan(SC_common_name):
                    SC_common_names.append(gene)
                else: 
                    print('float but not nan uh oh!')
            else: 
                SC_common_names.append(SC_common_name_lookup[gene])
        except KeyError: 
            SC_common_names.append(gene)
    
    return SC_common_names

def SC_orf_lookup_by_name(name_list):
    #Input is a list of common names, output is a list of orfs
    #would be nice if this worked with a single string as well as a list 
    #e.g. 'ERV14' as well as ['ERV14']
    
    SC_orfs_lookup, SC_common_name_lookup, SC_features_lookup = read_SGD_features()
    sc_genenames = []
    
    for name in name_list: 
        try:  
            sc_genename = SC_orfs_lookup[name]
        except KeyError: 
            print("S.Cer orf for " + name + "not found")
            sc_genename = name
        
        sc_genenames.append(sc_genename)
    
    return sc_genenames

def parse_pathway_list(pathway_list_str):
    if len(pathway_list_str)>0:
        if pathway_list_str[0]=='[':
            #parse list: 
            pathway_list = pathway_list_str.strip("[]'").split(',')
            if len(pathway_list)>1:
                print('Warning: more than one pathway for a gene')
            if pathway_list[0]=='':
                print('No_pathway in pathway_list_str, returning None')
                return None
            else: 
                return pathway_list[0]
        else: 
            print('pathway_list_str not a list, returning None')
            return None
    else:
        print('pathway_list_str is empty string, returning None')
        return None

def seq_record_fasta_printout(seq_records, f_out, gene_full_set, seqs_to_get, proteome_source, spec_orig_genome):
    
    for seq_record in seq_records:
        gene_full_check = gene_id_retrieve(proteome_source, seq_record)
        if (gene_full_check in gene_full_set):
            
            if proteome_source == 'shen': 
                (y1000_id, og, maxval, diff_top2_val) = seqs_to_get[gene_full_check]
                gene_full_shen = gene_full_check
            else: 
                (y1000_id, og, maxval, diff_top2_val, gene_full_shen) = seqs_to_get[gene_full_check]
                
            protein_seq = seq_record.seq
            f_out.write('>' + spec_orig_genome + '__' + og + '__' + y1000_id + 
                        ' source=' + proteome_source +
                        ' gene_full=' + gene_full_check +
                        ' gene_full_shen=' + gene_full_shen + 
                        ' L=' + str(len(protein_seq)) + 
                        ' sim_score_vs_shen={:0.1f}'.format(maxval) + 
                        ' sim_score_vs_shen_diff={:0.1f}'.format(diff_top2_val) + '\n')
            f_out.write(str(protein_seq) + '\n')  #I wonder why some of the bases were in lower case

def fasta_extract_shen(f_out, protein_dir, spec_orig_genome, y1000plus_dir, og_out_data, spec_id):  
    
    source = 'shen'
    sim_score_vs_shen = 'NA'
    sim_score_vs_shen_diff = 'NA'
    
    protein_fname = protein_dir + spec_orig_genome + '.max.pep'
    seq_records = SeqIO.parse(protein_fname, "fasta")

    spec_lookup_fname = y1000plus_dir + os.sep + os.path.normpath('y1000plus_tools_data/y1000plus/id_lookups/' + spec_orig_genome + '.csv')
    spec_lookup = pd.read_csv(spec_lookup_fname, index_col=0)

#             if spec_orig_genome == 'saccharomyces_cerevisiae':
#                 #for S. cerevisiae the orf name is the 'gene_full'
#                 gene_full_y1000_id_lookup = dict(zip(spec_lookup['y1000_id'], spec_lookup.index))
#             else: 
    gene_full_y1000_id_lookup = dict(zip(spec_lookup['y1000_id'], spec_lookup['gene_full']))

    #for each orthogroup:
    #extract protein seq if it comes from the right species.
    #gene_full: (y1000_id, og)
    seqs_to_get = {}
    gene_full_set = []
    for og, (N_genes, pct_present, og_genes_out) in og_out_data.items():
        for y1000_id in og_genes_out: 
            spec_id_check = int(y1000_id.split('_')[0])
            if spec_id_check==spec_id:
                gene_full = gene_full_y1000_id_lookup[y1000_id]
                seqs_to_get[gene_full] = (y1000_id, og)
                gene_full_set.append(gene_full)

    #gene_full_sets[spec_orig_genome] = gene_full_set
    
    #Go back through the protein fasta and if the gene is in gene_full_set, print the fasta line
    for seq_record in seq_records:
        #gene_full = 'augustus_masked-Deha2C-processed-gene-4.36'
        if spec_orig_genome == 'saccharomyces_cerevisiae':
            gene_full_check = seq_record.description.split()[0] #SC specific
        elif spec_orig_genome == 'candida_albicans': 
            gene_full_check = seq_record.description
        else: 
            gene_full_check = seq_record.description.split()[1].split('=')[1]
        #print(gene_full)
        if (gene_full_check in gene_full_set):
            ##find which y1000_id was matched
            #y1000_rlookup = genes_lookup['gene_full'] == gene_full
    #         for gene, tf in y1000_rlookup.items(): 
    #             if tf:
    #                 y1000_id=gene
    #         gene_id = genes_lookup.loc[y1000_id, 'gene_id']
            (y1000_id, og) = seqs_to_get[gene_full_check]
            protein_seq = seq_record.seq
            f_out.write('>' + spec_orig_genome + '__' + og + '__' + y1000_id + 
                        ' source=' + source +
                        ' gene_full=' + gene_full_check +
                        ' gene_full_shen=' + gene_full_check + 
                        ' L=' + str(len(protein_seq)) + 
                        ' sim_score_vs_shen=' + sim_score_vs_shen + 
                        ' sim_score_vs_shen_diff=' + sim_score_vs_shen_diff + '\n')
            f_out.write(str(protein_seq) + '\n')  #I wonder why some of the bases were in lower case
            
def load_model_protein_dict(spec_abbrev):
    #Load peptide sequences for model species, make dictionary from gene id to peptide sequence
    
    prot_fnames = {'Scer': genomes_dir +os.sep +  os.path.normpath('saccharomyces_cerevisiae\S288C_reference_genome_R64-2-1_20150113\orf_trans_all_R64-2-1_20150113.fasta'),
                   'Calb': genomes_dir +os.sep +  os.path.normpath('candida_albicans\C_albicans_SC5314_A22_current_default_protein.fasta'),
                   'Spom': genomes_dir +os.sep +  os.path.normpath('schizosaccharomyces_pombe\peptide_20240918.fa') #Note: previously may have used a version from 20221012
                  }

    prot_fname = prot_fnames[spec_abbrev]
    proteins = SeqIO.parse(prot_fname, 'fasta')
    proteins_seqs = {}
    for record in proteins: 
        if spec_abbrev == 'Spom':
            gene_id = '.'.join(record.id.split(':')[0].split('.')[0:2])  #  Go from e.g. 'SPAC1002.01.1:pep' to 'SPAC1002.01'
        else: 
            gene_id = record.id
        proteins_seqs[gene_id] = str(record.seq)
    
    return proteins_seqs

def load_model_gene_id_2_y1000_id():
    #ID Mapping for model species
    #Lookup from gene_id to y1000_id 
    gene_id_2_y1000_id = {}

    #Load S.cer lookup table: 
    scer_lookup_fname = y1000plus_dir + os.sep + os.path.normpath('y1000plus_tools_data/y1000plus/id_lookups/saccharomyces_cerevisiae.csv')
    scer_lookup = pd.read_csv(scer_lookup_fname, index_col=0)
    gene_id_2_y1000_id['Scer'] = dict(zip(scer_lookup.index,scer_lookup['y1000_id']))

    #Load C.alb lookup table
    calb_lookup_fname = y1000plus_dir + os.sep + os.path.normpath('y1000plus_tools_data/y1000plus/id_lookups/candida_albicans.csv')
    calb_lookup = pd.read_csv(calb_lookup_fname, index_col=0)
    gene_id_2_y1000_id['Calb'] = dict(zip(calb_lookup.index,calb_lookup['y1000_id']))

    #Spom is difficult because the y1000 id map goes from gene names rather than systematic ids and many of the gene names are synonyms rather than the official names
    #Load S.pom lookup table: 
    spom_lookup_fname = y1000plus_dir + os.sep + os.path.normpath('y1000plus_tools_data/y1000plus/id_lookups/schizosaccharomyces_pombe.csv')
    spom_lookup = pd.read_csv(spom_lookup_fname, index_col=0)
    spom_gene_ids = pd.read_table(genomes_dir + os.sep + os.path.normpath('schizosaccharomyces_pombe/gene_IDs_names_products.tsv'), header=None)
    gene_full_2_y1000_id = dict(zip(spom_lookup['gene_full'],spom_lookup['y1000_id']))
    spom_name_2_systematic_id = dict(zip(spom_gene_ids[2],spom_gene_ids[0]))

    #map synonyms
    syn_map_2_name = {}
    syn_map_2_systematic_id = {}
    for syn_ids, name, systematic_id in zip(spom_gene_ids[7], spom_gene_ids[2], spom_gene_ids[0]): 
        if not(pd.isna(syn_ids)): 
            syn_ids_split = syn_ids.split(',')

            for syn_id in syn_ids_split: 
                syn_map_2_name[syn_id] = name
                syn_map_2_systematic_id[syn_id] = systematic_id

    gene_id_2_y1000_id['Spom'] = {}
    for gene_full, y1000_id in gene_full_2_y1000_id.items():
        if gene_full in set(spom_gene_ids[0]): #gene_full is a systematic id
            gene_id_2_y1000_id['Spom'][gene_full] = y1000_id
        elif gene_full in spom_name_2_systematic_id.keys(): 
            systematic_id = spom_name_2_systematic_id[gene_full]
            gene_id_2_y1000_id['Spom'][systematic_id] = y1000_id
        elif gene_full in syn_map_2_systematic_id.keys(): 
            #print('Exception gene for pombe: ' + gene_full)
            systematic_id = syn_map_2_systematic_id[gene_full]
            if (spom_name_2_systematic_id[syn_map_2_name[gene_full]]!=systematic_id): 
                print('gene_full not mapping to name that maps to systematic id ' + gene_full + ' systematic_id=' + systematic_id)
            gene_id_2_y1000_id['Spom'][systematic_id] = y1000_id
        else: 
            print('Spom: No official name, systematic_id or synonym for ' + gene_full + ' y1000_id=' + y1000_id)

    return gene_id_2_y1000_id

def load_model_swissprot_id_2_gene_id():
    #Lookup from swissprot_id to gene_id
    #would prefer to use a more official source for these for S.cer and S.pom

    swissprot_id_2_gene_id= {}
    scer_swissprot_id_2_gene_id_df = pd.read_table(base_dir + os.sep + os.path.normpath('msas/structural/Scer_protein_names.tsv'))
    swissprot_id_2_gene_id['Scer'] = dict(zip(scer_swissprot_id_2_gene_id_df['Swiss-Prot'],scer_swissprot_id_2_gene_id_df['OLN']))
    
    cgd_id_2_swissprot_id_df = pd.read_table(genomes_dir + os.sep + 'candida_albicans' + os.sep + 'gp2protein_C_albicans_SC5314', header=None)
    cgd_id_2_swissprot_id_df['uniprot']=[uniprot_id.split(':')[1] for uniprot_id in cgd_id_2_swissprot_id_df[1]]
    cgd_id_2_swissprot_id_df['cgd']=[cgd_id.split(':')[1] for cgd_id in cgd_id_2_swissprot_id_df[0]] 
    #cgd_id_2_swissprot_id = dict(zip(cgd_id_2_swissprot_id_df['cgd'],cgd_id_2_swissprot_id_df['uniprot']))

    calb_table = pd.read_table(genomes_dir + os.sep + 'candida_albicans' + os.sep + 'C_albicans_SC5314_version_A22-s07-m01-r164_chromosomal_feature.tab', skiprows=8, header=None)
    #gene_id_2_cgd_id = dict(zip(calb_table[0], calb_table[8]))
    calb_gene_id_2_swissprot_id_df = pd.merge(calb_table.loc[:,[0,8]],  cgd_id_2_swissprot_id_df.loc[:,['uniprot','cgd']], how='left', left_on=8, right_on='cgd')
    calb_gene_id_2_swissprot_id_df_uniprot = calb_gene_id_2_swissprot_id_df.loc[[not(item) for item in list(calb_gene_id_2_swissprot_id_df['uniprot'].isna())],:]

    swissprot_id_2_gene_id['Calb'] = dict(zip(calb_gene_id_2_swissprot_id_df_uniprot['uniprot'], calb_gene_id_2_swissprot_id_df_uniprot[0]))

    spom_swissprot_id_2_gene_id_df = pd.read_table(genomes_dir + os.sep + os.path.normpath('schizosaccharomyces_pombe/PomBase2UniProt.tsv'), header=None)
    swissprot_id_2_gene_id['Spom'] = dict(zip(spom_swissprot_id_2_gene_id_df[1],spom_swissprot_id_2_gene_id_df[0]))
    
    return swissprot_id_2_gene_id

def load_model_gene_id_2_swissprot_id():
    #Lookup from gene_id to swissprot_id (uniprot). 
    #would prefer to use a more official source for these for S.cer and S.pom

    gene_id_2_swissprot_id= {}
    scer_gene_id_2_swissprot_id_df = pd.read_table(base_dir + os.sep + os.path.normpath('msas/structural/Scer_protein_names.tsv'))
    gene_id_2_swissprot_id['Scer'] = dict(zip(scer_gene_id_2_swissprot_id_df['OLN'], scer_gene_id_2_swissprot_id_df['Swiss-Prot']))

    
    cgd_id_2_swissprot_id_df = pd.read_table(genomes_dir + os.sep + 'candida_albicans' + os.sep + 'gp2protein_C_albicans_SC5314', header=None)
    cgd_id_2_swissprot_id_df['uniprot']=[uniprot_id.split(':')[1] for uniprot_id in cgd_id_2_swissprot_id_df[1]]
    cgd_id_2_swissprot_id_df['cgd']=[cgd_id.split(':')[1] for cgd_id in cgd_id_2_swissprot_id_df[0]] 
    #cgd_id_2_swissprot_id = dict(zip(cgd_id_2_swissprot_id_df['cgd'],cgd_id_2_swissprot_id_df['uniprot']))

    calb_table = pd.read_table(genomes_dir + os.sep + 'candida_albicans' + os.sep + 'C_albicans_SC5314_version_A22-s07-m01-r164_chromosomal_feature.tab', skiprows=8, header=None)
    #gene_id_2_cgd_id = dict(zip(calb_table[0], calb_table[8]))
    calb_gene_id_2_swissprot_id_df = pd.merge(calb_table.loc[:,[0,8]],  cgd_id_2_swissprot_id_df.loc[:,['uniprot','cgd']], how='left', left_on=8, right_on='cgd')
    gene_id_2_swissprot_id['Calb'] = dict(zip(calb_gene_id_2_swissprot_id_df[0], calb_gene_id_2_swissprot_id_df['uniprot']))

    spom_gene_id_2_swissprot_id_df = pd.read_table(genomes_dir + os.sep + os.path.normpath('schizosaccharomyces_pombe/PomBase2UniProt.tsv'), header=None)
    gene_id_2_swissprot_id['Spom'] = dict(zip(spom_gene_id_2_swissprot_id_df[0], spom_gene_id_2_swissprot_id_df[1]))
    
    return gene_id_2_swissprot_id

def seq_squeeze(seq_in, all_outputs=True):
    #input: A string containing a sequence with dashes from a multiple sequence alignment
    #Output: 
    # seq_out - a string with a squeezed sequence with dashes removed
    # msa2input: a dictionary mapping indices from MSA to indices for the squeezed sequence (0 based)
    # input2msa: a dictionary mapping indices from input sequence to the MSA (0 based)
    # pair_mapping - a mapping containing consecutive integers that link the index from the squeezed sequence to the corresponding indices of the msa.
    # 
    
    seq_out = ''
    for res in seq_in:
        if res!= '-':
            seq_out = seq_out + res
    
    if all_outputs: 
        msa2input = {}
        input2msa = {}
        seq_ind = 0
        #msa_ind_old = 0
        pair_mapping = []
        for msa_ind, res in enumerate(seq_in):
            if res != '-':
                msa2input[msa_ind] = seq_ind
                input2msa[seq_ind] = msa_ind
                if seq_ind>0:
                    seq_pair = (seq_ind-1,seq_ind)
                    msa_pair = (msa_ind_old,msa_ind)
                    pair_mapping.append((seq_pair, msa_pair))

                seq_ind = seq_ind + 1
                msa_ind_old = msa_ind
                #print(msa_ind_old)
                #print(seq_ind)

            else: 
                msa2input[msa_ind] = None

    
        return seq_out, msa2input, input2msa, pair_mapping
    
    else:
        return seq_out

def load_model_og_lookup():
    #Make dictionary to look up og from gene_id for model species from og_metadatafile
    
    og_metadata_fname = base_dir + os.sep + os.path.normpath('selected_proteins/og_metadata.json')

    with open(og_metadata_fname, 'r') as f:
        og_metadata = json.load(f) 
    
    #Make look up table from og_metadata
    abbrev_fields = [('Calb', 'calb_genes'),
                     ('Scer', 'sc_genes'), 
                     ('Spom', 'spom_genes')]

    model_og_lookup = {}

    for (spec_abbrev, metadata_field) in abbrev_fields: 
        model_og_lookup_spec = {}
        for og, metadata in og_metadata.items(): 
            model_genes = metadata[metadata_field]
            for gene in model_genes:
                if gene != 'NONE': 
                    model_og_lookup_spec[gene] = og
        model_og_lookup[spec_abbrev] = model_og_lookup_spec
    
    return model_og_lookup

def calMean_dN_dS_yn00(paml_gene_dn_ds_file, output_file, method="median dN_dS, dN, dS", min_ds = 0.005, max_ds = 3.0, max_dn_ds = 50):
    #Parse YN00 files 
    #Based on https://github.com/SysBioChalmers/Multi_scale_evolution/blob/main/evolution_analysis/code/gene_dn_ds_paml/result_parse_yn00_update.py
    
    result_file = open(paml_gene_dn_ds_file).readlines()
    index0 = [i for i, x in enumerate(result_file) if "seq. seq." in x]
    if len(index0)==0: 
        warnings.warn('No dN dS table in results file for ' + paml_gene_dn_ds_file)
        if method in ['mean', 'median', 'max']:
            return np.nan
        elif method=='median dN_dS, dN, dS':
            return (np.nan, np.nan, np.nan)
    else: 
    
        index1 = [i for i, x in enumerate(result_file) if "(C) LWL85, LPB93 & LWLm methods" in x]
        dn_ds = result_file[index0[0]:index1[0]]
        # further remove the line with only "\n"
        dn_ds0 = [x for x in dn_ds if x != "\n"]
        # save the dn_ds
        dn_ds_df = Export_dn_ds_yn00(dnds1=dn_ds0, output_file=output_file)

        #Filter dS to be between 0.005 and less than 3
        dn_ds_filt1 = dn_ds_df[(dn_ds_df['dS_new']>min_ds)&(dn_ds_df['dS_new']<max_ds)]
        if len(dn_ds_filt1)==0: 
            warnings.warn('No dS values pass filters (min_ds = ' + str(min_ds) + ', max_ds = ' + str(max_ds) + ') for ' + paml_gene_dn_ds_file)

        #Filter Omega greater than 50 
        #Original code had this comment:
        #       here should we filter out omega > 50
        #       here it may be not reasonable to filter out omega > 50. In fact, these value should be zero
        # I presume the DS filter gets rid of those, i will raise a warning if this filter is used
        dn_ds_filt2 = dn_ds_filt1[dn_ds_filt1['omega_new']<max_dn_ds]
        if len(dn_ds_filt2)<len(dn_ds_filt1): 
            warnings.warn('Omega greater than ' + str(max_dn_ds) + ' in ' + paml_gene_dn_ds_file)

        ds = dn_ds_filt2['dS_new']
        dn = dn_ds_filt2['dN_new']
        dn_ds = dn_ds_filt2['omega_new']


        if method == "mean":
            average_dn_ds = np.mean(dn_ds[np.isfinite(dn_ds)])
            return average_dn_ds
        elif method == "median":
            median_dn_ds = np.median(dn_ds[np.isfinite(dn_ds)])
            return median_dn_ds
        elif method == "max":
            max_dn_ds = np.nanmax(dn_ds)
            return max_dn_ds
        elif method == "median dN_dS, dN, dS":
            return (np.median(dn_ds[np.isfinite(dn_ds)]),np.median(dn[np.isfinite(dn)]),np.median(ds[np.isfinite(ds)]) )
            
def Export_dn_ds_yn00(dnds1, output_file):
    new_dnds = []
    colname1 = None
    for i, x in enumerate(dnds1):
        if i == 0:
            colname = x.split("  ")
            colname0 = [n.strip(" ").strip("\n") for n in colname if n != ""]
            c1 = ["seq1", "seq2"]
            colname1 = c1 + colname0[1:]
        else:
            x1 = x.split(" ")
            x1 = [n.strip(" ").strip("\n") for n in x1 if n != ""]
            n1 = " ".join(x1[10:])
            n2 = " ".join(x1[7:10])
            x1_new = x1[0:7] + [n2] + [n1]
            x2 = "@".join(x1_new)
            new_dnds.append(x2)
    df = pd.DataFrame({"ID": new_dnds})
    df1 = df['ID'].str.split('@', expand=True)
    df1.columns = colname1

    dS0 = df1["dS +- SE"].tolist()
    dS1 = [float(x.split("+-")[0]) for x in dS0]
    dN0 = df1["dN +- SE"].tolist()
    dN1 = [float(x.split("+-")[0]) for x in dN0]
    omega0 = df1["omega"].tolist()
    omega1 = [float(x) for x in omega0]
    

    df1["dS_new"] = dS1
    df1["dN_new"] = dN1
    df1["omega_new"] = omega1

    df1.to_csv(output_file)
    return df1

def load_m0_data(m0_base, branch_table_out=False):
    #Extracts m0 data from a folder containing folders with output from m0 calculations. 
    #m0_base = base_dir + os.sep + os.path.normpath('selection_calculations/m0') + os.sep
    #branch_table_out = True

    og_list = []
    og_ref_list = []
    m0_data_out = {}

    if branch_table_out: 
        branch_tables = {}

    for og_ref in next(os.walk(m0_base))[1]:    #next(os.walk(m0_base))[0] gets just directories - https://stackoverflow.com/questions/141291/how-to-list-only-top-level-directories-in-python
        ##for og_ref in tm_align_post_trim_filter_list:
        #og_ref = 'OG1299_REF_Scer_AF-P00549-F1-model_v2'
        og = og_ref.split('_')[0]
        if og[0:2] =='OG':
            og_ref_list.append(og_ref)
            og_list.append(og)
            #print(og_ref)
            m0_dir = m0_base + og_ref + os.sep
            output_file = m0_dir + 'm0.csv'
            paml_gene_dn_ds_file = m0_dir +  'm0.out'

            columns_base = ['og', 'tree_length', 'kappa', 'dN_dS_struct', 'tree_length_dN', 'tree_length_dS', 'convergence_issue', 'dN_dS_found']
            tree_length = None
            kappa = None
            dN_dS = None
            tree_length_dN = None
            tree_length_dS=None
            convergence_issue = None

            #Branch table
            branch_table = None

            #Branch_table summary quantities: 
            branch_table_summary_columns = ['dN_mean', 'dN_std', 'dS_mean', 'dS_std']
            dN_mean = None
            dN_std = None
            dS_mean = None
            dS_std = None

            with open(paml_gene_dn_ds_file, 'r') as m0_out_data: 
                for line in m0_out_data:
                    if 'TREE #' in line: 
                        line = next(m0_out_data)
                        convergence_issue = False
                        if 'check convergence..' in line:
                            print('convergence issue in ' + og_ref)
                            convergence_issue = True              
                        break

                for line in m0_out_data:
                    if 'tree length =' in line: 
                        tree_length = float(line.split('=')[1].strip())
                        break

                dN_dS_found = True
                #This assumes that if there is no tree lenth found, then there is no further calculation
                if tree_length==None:
                    dN_dS_found = False
                    print('dN_dS not calculated for' + og_ref)


                if dN_dS_found: 

                    for line in m0_out_data:
                        if 'kappa' in line: 
                            kappa = float(line.split('=')[1].strip())
                            break

                    for line in m0_out_data:
                        if 'omega' in line: 
                            dN_dS = float(line.split('=')[1].strip())
                            break

                    if branch_table_out: 
                        for line in m0_out_data: 
                            if 'dN & dS for each branch' in line: 
                                next(m0_out_data)
                                line = next(m0_out_data)
                                branch_table_headers = line.split()
                                next(m0_out_data)
                                line = next(m0_out_data)
                                linesp = line.split()
                                branch_table_data = []
                                while len(linesp)==len(branch_table_headers): 
                                    branch_table_data.append(linesp)
                                    line = next(m0_out_data)
                                    linesp=line.split()
                                break
                        branch_table = pd.DataFrame(branch_table_data, columns = branch_table_headers)
                        #Make all columns except 'branch' numeric columns
                        branch_table_headers_numeric = branch_table_headers.copy()
                        branch_table_headers_numeric.remove('branch')
                        for col in branch_table_headers_numeric: 
                            branch_table[col] = pd.to_numeric(branch_table[col])
                        dN_mean=np.mean(branch_table['dN'])
                        dN_std=np.std(branch_table['dN'])
                        dS_mean=np.mean(branch_table['dS'])
                        dS_std=np.std(branch_table['dS'])

                    for line in m0_out_data:
                        if 'tree length for dN:' in line: 
                            tree_length_dN = float(line.split(':')[1].strip())
                            break

                    for line in m0_out_data:
                        if 'tree length for dS:' in line: 
                            tree_length_dS = float(line.split(':')[1].strip())
                            break

                if branch_table_out: 
                    m0_data_out[og_ref] = (og, tree_length, kappa, dN_dS, tree_length_dN, tree_length_dS, convergence_issue, dN_dS_found, dN_mean, dN_std, dS_mean, dS_std)
                    branch_tables[og_ref] = branch_table
                else: 
                    m0_data_out[og_ref] = (og, tree_length, kappa, dN_dS, tree_length_dN, tree_length_dS, convergence_issue, dN_dS_found)
        else: 
            print('Folder ' + og_ref + ' does not contain dnds calculation')

    if branch_table_out: 

        m0_data_df = pd.DataFrame.from_dict(m0_data_out, orient='index', columns = columns_base + branch_table_summary_columns)
        return m0_data_df, branch_tables

    else: 
        m0_data_df = pd.DataFrame.from_dict(m0_data_out, orient='index', columns = columns_base)
        return m0_data_df

def extract_beb_values(bs_rst_fname):
    #Extract Bayes Empirical Bayes Values for each residue in a Branch/Site test calculation
    site_classes = {1:'class_0', 
                2:'class_1', 
                3:'class_2a',
                4:'class_2b'
               }

    neb_data = {}

    with open(bs_rst_fname, 'r') as f_in: 
        #result_file = open(bs_rst_file).readlines()

        for line in f_in:
            if 'TREE #' in line: 
                tree_no = int(line.split()[2])
                neb_data_tree = {}

                beb_start = False
                while not(beb_start):
                    line = next(f_in)
                    if "Bayes Empirical Bayes" in line:
                        beb_start=True

                next(f_in)
                next(f_in)

                line = next(f_in)
                line_sp = line.split()
                while len(line_sp)==8:        
                    (ind, aa, class_0, class_1, class_2a, class_2b, paren, class_ind_paren) = line_sp
                    neb_data_tree[int(ind)] = (aa, float(class_0), float(class_1), float(class_2a), float(class_2b), float(class_2a)+float(class_2b), site_classes[int(class_ind_paren.split(')')[0])])
                    line = next(f_in)
                    line_sp = line.split()

                neb_data[tree_no] = pd.DataFrame.from_dict(neb_data_tree, orient='index', columns=['ref_AA','class_0','class_1','class_2a','class_2b','class_2','pred_class'])

    return neb_data

def load_site_model_data(site_model_base, branch_table_out=False):
    #Extracts output data from codeml from a folder containing folders with output from calculations using NSsites Models 0 (m0), 1 (NearlyNeutral), 3(PositiveSelection) calculations. 
    #site_model_base = base_dir + os.sep + os.path.normpath('selection_calculations/site_model') + os.sep
    #branch_table_out = True

    og_list = []
    og_ref_list = []

    tests = ['m0','NearlyNeutral','PositiveSelection']

    test_nclasses = {'NearlyNeutral': 2,
                    'PositiveSelection': 3
                    }

    test_line_text = {'m0': 'NSsites Model 0: one-ratio', 
                    'NearlyNeutral': 'NSsites Model 1: NearlyNeutral (2 categories)' ,
                    'PositiveSelection': 'NSsites Model 2: PositiveSelection (3 categories)'}


    columns_base = {'m0': ['og', 'mle', 'tree_length', 'kappa', 'dN_dS', 'tree_length_dN', 'tree_length_dS', 'convergence_issue', 'dN_dS_found'], 
                    'NearlyNeutral': ['og','mle','tree_length', 'kappa', 'p_1', 'dN_dS_1', 'p_2', 'dN_dS_2', 'convergence_issue', 'dN_dS_found'],
                    'PositiveSelection': ['og','mle','tree_length', 'kappa', 'p_1', 'dN_dS_1', 'p_2', 'dN_dS_2', 'p_3', 'dN_dS_3', 'convergence_issue', 'dN_dS_found']
                    }

    site_model_data_out = {test:{} for test in tests}
    site_model_data_dfs = {test:{} for test in tests}

    if branch_table_out: 
        branch_tables = {test:{} for test in tests}

    for og_ref in next(os.walk(site_model_base))[1]:    #next(os.walk(m0_base))[0] gets just directories - https://stackoverflow.com/questions/141291/how-to-list-only-top-level-directories-in-python
        #         ##for og_ref in tm_align_post_trim_filter_list:

        #og_ref = 'OG4555_REF_Scer_AF-P32799-F1-model_v2'  #example
        og = og_ref.split('_')[0]
        if og[0:2] =='OG':
            og_ref_list.append(og_ref)
            og_list.append(og)
            #print(og_ref)
            site_model_dir = site_model_base + og_ref + os.sep
            output_files = {test:{site_model_dir + 'site_model_' + test + '.csv'} for test in tests},
            paml_gene_dn_ds_file = site_model_dir +  'sm01278.out' #Note: filename is 01278, but most of the site model runs were just with models 0,1,2 because model 8 especially kept going into an infinite loop

            for test in tests:
            
                #test = 'PositiveSelection'
                site_model_data_out_test = site_model_data_out[test]
                if branch_table_out: 
                    branch_tables_test = branch_tables[test]


                columns_base_test = columns_base[test]

                data_out_test = {column: None for column in columns_base_test}
                data_out_test['og']=og

                # tree_length = None
                # kappa = None
                # dN_dS = None
                # tree_length_dN = None
                # tree_length_dS=None
                # convergence_issue = None

                #Branch table
                branch_table = None

                #Branch_table summary quantities: 
                branch_table_summary_columns = ['dN_mean', 'dN_std', 'dS_mean', 'dS_std']
                dN_mean = None
                dN_std = None
                dS_mean = None
                dS_std = None

                with open(paml_gene_dn_ds_file, 'r') as out_data: 
                    for line in out_data: 
                        if test_line_text[test] in line: 
                            break
                            #line=next(out_data)
                    
                    for line in out_data:
                        if 'TREE #' in line: 
                            line = next(out_data)
                            data_out_test['convergence_issue'] = False
                            if 'check convergence..' in line:
                                print('convergence issue in ' + og_ref +' for ' + test + ' test.')
                                data_out_test['convergence_issue'] = True              
                                line = next(out_data)
                            #get maximum likelihood estimate
                            data_out_test['mle'] = float(line.split(':')[3].split()[0])
                            break

                    for line in out_data:
                        if 'tree length =' in line: 
                            data_out_test['tree_length'] = float(line.split('=')[1].strip())
                            break

                    data_out_test['dN_dS_found'] = True
                    #This assumes that if there is no tree length found, or if the mle is nan, then there is no further calculation
                    if data_out_test['tree_length']==None:
                        data_out_test['dN_dS_found']= False
                        print('dN_dS not calculated for ' + og_ref + ' test ' + test + ': Zero Tree length')
                    elif np.isnan(data_out_test['mle']):
                        data_out_test['dN_dS_found']= False
                        print('dN_dS not calculated for ' + og_ref + ' test ' + test + ': mle is nan')

                    if data_out_test['dN_dS_found']: 

                        for line in out_data:
                            if 'kappa' in line: 
                                data_out_test['kappa'] = float(line.split('=')[1].strip())
                                break

                        if test=='m0':
                            for line in out_data:
                                if 'omega' in line: 
                                    data_out_test['dN_dS'] = float(line.split('=')[1].strip())
                                    break
                        elif test in ['NearlyNeutral','PositiveSelection']:
                            n_classes = test_nclasses[test]
                            for line in out_data: 
                                if 'p:' in line:  
                                    for jj in range(1,n_classes + 1):
                                        data_out_test['p_' + str(jj)] = float(line.split()[jj].strip())
                                    break    

                            for line in out_data: 
                                if 'w:' in line:  
                                    for jj in range(1,n_classes + 1):
                                        data_out_test['dN_dS_' + str(jj)] = float(line.split()[jj].strip())
                                    break
                    

                        if branch_table_out: 
                            for line in out_data: 
                                if 'dN & dS for each branch' in line: 
                                    next(out_data)
                                    line = next(out_data)
                                    branch_table_headers = line.split()
                                    next(out_data)
                                    line = next(out_data)
                                    linesp = line.split()
                                    branch_table_data = []
                                    while len(linesp)==len(branch_table_headers): 
                                        branch_table_data.append(linesp)
                                        line = next(out_data)
                                        linesp=line.split()
                                    break
                            branch_table = pd.DataFrame(branch_table_data, columns = branch_table_headers)
                            #Make all columns except 'branch' numeric columns
                            branch_table_headers_numeric = branch_table_headers.copy()
                            branch_table_headers_numeric.remove('branch')
                            for col in branch_table_headers_numeric: 
                                branch_table[col] = pd.to_numeric(branch_table[col])
                            dN_mean=np.mean(branch_table['dN'])
                            dN_std=np.std(branch_table['dN'])
                            dS_mean=np.mean(branch_table['dS'])
                            dS_std=np.std(branch_table['dS'])

                        if test=='m0':
                            for line in out_data:
                                if 'tree length for dN:' in line: 
                                    data_out_test['tree_length_dN'] = float(line.split(':')[1].strip())
                                    break

                            for line in out_data:
                                if 'tree length for dS:' in line: 
                                    data_out_test['tree_length_dS'] = float(line.split(':')[1].strip())
                                    break

                    if branch_table_out: 
                        site_model_data_out_test[og_ref] = [value for key, value in data_out_test.items()] + [dN_mean, dN_std, dS_mean, dS_std]
                        #og, tree_length, kappa, dN_dS, tree_length_dN, tree_length_dS, convergence_issue, dN_dS_found, dN_mean, dN_std, dS_mean, dS_std)
                        branch_tables_test[og_ref] = branch_table
                    else: 
                        site_model_data_out_test[og_ref] = [value for key, value in data_out_test.items()]
                        #(og, tree_length, kappa, dN_dS, tree_length_dN, tree_length_dS, convergence_issue, dN_dS_found)
                site_model_data_out[test] = site_model_data_out_test
                if branch_table_out: 
                    branch_tables[test] = branch_tables_test
        else: 
            print('Folder ' + og_ref + ' does not contain dnds calculation')

    if branch_table_out: 
        for test in tests: 
            site_model_data_dfs[test] = pd.DataFrame.from_dict(site_model_data_out[test], orient='index', columns = columns_base[test] + branch_table_summary_columns)
        return site_model_data_dfs, branch_tables

    else: 
        for test in tests: 
            site_model_data_dfs[test] = pd.DataFrame.from_dict(site_model_data_out[test], orient='index', columns = columns_base[test])
        return site_model_data_dfs

def parse_sitemodel_NEB_2(line):
    #example line: 1 A   0.00000 1.00000 ( 2)  1.000'
    linesp = line.split()
    output =  (int(linesp[0]), linesp[1], float(linesp[2]), float(linesp[3]),linesp[5].strip(')'),float(linesp[6]))
    return output

def parse_sitemodel_NEB_3(line):
    #example line: 1 A   0.00000 1.00000 0.00000 ( 2)  1.000  0.000
    linesp = line.split()
    if len(linesp)==9:  #If it passes the positive seleciton test, it looks like it calculates an extra value
        output =  (int(linesp[0]), linesp[1], float(linesp[2]), float(linesp[3]), float(linesp[4]), linesp[6].strip(')'),float(linesp[7]),float(linesp[8]))
    elif len(linesp)==8: 
        output =  (int(linesp[0]), linesp[1], float(linesp[2]), float(linesp[3]), float(linesp[4]), linesp[6].strip(')'),float(linesp[7]),np.nan)
    return output

def parse_sitemodel_BEB_3(line):
    #example line: 1 A   0.00000 1.00000 0.00000 ( 2)  1.000 +-  0.000
    linesp = line.split()
    output =  (int(linesp[0]), linesp[1], float(linesp[2]), float(linesp[3]), float(linesp[4]), linesp[6].strip(')'),float(linesp[7]), float(linesp[9]))
    return output

def extract_posterior_values_site_model(og_ref, site_model_dir, save_output=False):
    #output is three tables.  NEB table for 2 classes, NEB table for 2 classes, BEB table for 3 classes

    bs_rst_fname = site_model_dir + os.sep + og_ref + os.sep + 'rst'

    #Extract Bayes Empirical Bayes Values for each residue in a Branch/Site test calculation

    site_classes = {'NEB_2': {
                        'classes' : ['1','2'], 
                        'header' : 'Naive Empirical Bayes (NEB) probabilities for 2 classes & postmean_w',
                        'columns' : ['strict_trim_ind','residue_first_seq','prob_purifying','prob_neutral','class','dN_dS'],
                        'parse_function': parse_sitemodel_NEB_2
                        },
                    'NEB_3': {
                        'classes' : ['1','2','3'], 
                        'header' : 'Naive Empirical Bayes (NEB) probabilities for 3 classes & postmean_w',
                        'columns' : ['strict_trim_ind','residue_first_seq','prob_purifying','prob_neutral', 'prob_positive','class','dN_dS','P(dN_dS>1)'],
                        'parse_function': parse_sitemodel_NEB_3
                        },
                    'BEB_3': {
                        'classes' : ['1','2','3'], 
                        'header' : 'Bayes Empirical Bayes (BEB) probabilities for 3 classes (class) & postmean_w',
                        'columns' : ['strict_trim_ind','residue_first_seq','prob_purifying','prob_neutral', 'prob_positive','class','dN_dS','P(dN_dS>1)'],
                        'parse_function': parse_sitemodel_BEB_3
                        }
                    }

    posterior_data = {}


    #Naive Empirical Bayes (NEB) probabilities for 2 classes & postmean_w
    #Naive Empirical Bayes (NEB) probabilities for 3 classes & postmean_w
    #Bayes Empirical Bayes (BEB) probabilities for 3 classes (class) & postmean_w
    #   1 A   0.00000 1.00000 0.00000 ( 2)  1.000 +-  0.000

    with open(bs_rst_fname, 'r') as f_in: 
        #result_file = open(bs_rst_file).readlines()

        for site_class in ['NEB_2','NEB_3','BEB_3']:  #need to iterate in correct order given file format
            (classes,header,columns, parse_function) = list(site_classes[site_class].values())
            outputs = []
            for line in f_in:
                #print(line)
                if header in line: 
                    next(f_in)
                    next(f_in)

                    line = next(f_in)
                    line_sp = line.split()
                                
                    while len(line_sp)>2:     
                        output = parse_function(line)
                        outputs.append(output)
                        line = next(f_in)
                        line_sp = line.split()
                    break
            
            outputs_df = pd.DataFrame(outputs, columns=columns)
            outputs_df.set_index('strict_trim_ind', inplace=True)

            posterior_data[site_class] = outputs_df

            if save_output:
                outputs_df.to_csv(site_model_dir + os.sep + og_ref + os.sep + site_class + '.csv', na_rep='nan')

            #Save outputs_df data as csv


            #         while len(line_sp)==8:        
            #             (ind, aa, class_0, class_1, class_2a, class_2b, paren, class_ind_paren) = line_sp
            #             neb_data_tree[int(ind)] = (aa, float(class_0), float(class_1), float(class_2a), float(class_2b), float(class_2a)+float(class_2b), site_classes[int(class_ind_paren.split(')')[0])])
            #             line = next(f_in)
            #             line_sp = line.split()

            #         neb_data[tree_no] = pd.DataFrame.from_dict(neb_data_tree, orient='index', columns=['ref_AA','class_0','class_1','class_2a','class_2b','class_2','pred_class'])

    return posterior_data

def identify_ref_residue(og_ref, residues_of_interest):
    #og_ref is the long name of the orthogroup you are interested in, e.g. 'OG1122_REF_Scer_AF-P13711-F1-model_v2'
    #residues of interest a list of tuples with the index and single letter AA code of the residue of interest from the strictly trimmed alignment

    ref = '_'.join(og_ref.split('_')[1:])

    #Load trimmed alignment
    aln_trim_fname = base_dir + os.sep + os.path.normpath('msas/structural/tm_align/trim_strict/' + og_ref + '.tm.fasta.clipkit')
    aln_trim = AlignIO.read(open(aln_trim_fname),'fasta')

    #Load clipkit log
    aln_trim_log_fname = aln_trim_fname + '.log'
    aln_trim_log = pd.read_table(aln_trim_log_fname, sep=' ', header=None, index_col=0)
    aln_trim_log.rename(columns={1:'trim',2:'note',3:'value'}, inplace=True)

    #Map from trimmed index to original index
    #This is a 1-based to 1-based map
    trim_ind_2_orig_ind = {}
    aln_trim_log_kept = aln_trim_log[aln_trim_log['trim']=='keep']
    for trim_ind, orig_ind in enumerate(aln_trim_log_kept.index): 
        trim_ind_2_orig_ind[trim_ind+1] = orig_ind

    aln_fname = base_dir + os.sep + os.path.normpath('msas/structural/tm_align/fasta_renamed/' + og_ref + '.tm.fasta')

    #Get species of first sequence (paml_spec used as the reference for the paml coordinates)
    #Also get index of Reference sequence
    ref_ind=None
    aln = AlignIO.read(open(aln_fname),'fasta')
    for (jj, record) in enumerate(aln): 
        if jj==0: 
            paml_spec = species_from_fasta_id(record.id)[0]
        #print(record.id.split('.')[0])
        if ref == record.id.split('.')[0]:
            ref_ind = jj

    #Get map for reference sequence from msa index
    ref_seq, msa2ref, ref2msa, pair_mapping = seq_squeeze(str(aln[ref_ind,:].seq))
    #msa2ref links the 0 based index of the alignment to the 0 based index of the reference sequence

    #after iteration need to read alignment in again
    aln = AlignIO.read(open(aln_fname),'fasta')


    ref_res_of_int = []
    aln_res_of_int = []

    for (res,aa) in residues_of_interest: 

        #aa from top line is 1-based index
        #aln_trim is 0-based
        if aa!=aln_trim[0,res-1]: 
            #Check if it is a CTG clade issue
            CTG_issue = False
            print(paml_spec)
            if paml_spec in ctg_clade_dict.keys(): 
                #If not raise an Error
                clade= ctg_clade_dict[paml_spec]
                alt_res = ctg_clade_alt_res[clade]
                
                if (aa == 'L') & (aln_trim[0,res-1]==alt_res):
                    CTG_issue = True
                    print('paml used wrong genetic code to identify residue in first species in alignment (' + paml_spec + ') ' + aa + str(res) + ' from paml should be ' + aln_trim[0,res-1])
                            
            if not(CTG_issue): 
                raise ValueError('Residue of first line does not match output from paml ' + aa + str(res) + ' != ' + aln_trim[0,res-1])

        #Identify index in original alignment
        #trim_ind_2_orig_ind is 1 based, orig_res_ind is 0-based
        orig_res_ind = trim_ind_2_orig_ind[res]-1

        orig_res_of_int = aln[0,orig_res_ind]
        #print(aln[:,(orig_res_ind-3):(orig_res_ind+3)])
        #print(orig_res_of_int)
        if aa!=orig_res_of_int: 
            #Check if it is a CTG clade issue
            CTG_issue = False
            if paml_spec in ctg_clade_dict.keys(): 
                #If not raise an Error
                clade= ctg_clade_dict[paml_spec]
                alt_res = ctg_clade_alt_res[clade]
                
                if (aa == 'L') & (aln_trim[0,res-1]==alt_res):
                    CTG_issue = True
                    print('paml used wrong genetic code to identify residue in first species in alignment (' + paml_spec + ') ' + aa + str(res) + ' from paml should be ' + orig_res_of_int)
                            
            if not(CTG_issue): 
                raise ValueError('Residue of first line in full alignment does not match output from paml ' + aa + str(res) + ' != ' + orig_res_of_int)

        #orig_res_ind is 0-based, msa2ref is 0 based, want ref_res_ind should be 1 based to match standard nomenclature
        if msa2ref[orig_res_ind] == None: 
            for ((ref_low,ref_high),(msa_low,msa_high)) in pair_mapping: 
                if (orig_res_ind > msa_low) & (orig_res_ind < msa_high): 
                    (ref_low_selected, ref_high_selected) = (ref_low, ref_high)
                    (msa_low_selected, msa_high_selected) = (msa_low, msa_high)
                    break
        
            ref_res_of_int.append('Gap between ' + aln[ref_ind,msa_low_selected] + str(ref_low_selected + 1) +' and ' + aln[ref_ind,msa_high_selected] + str(ref_high_selected + 1))
        
        else: 
            ref_res_ind = msa2ref[orig_res_ind] + 1
            ref_res_of_int.append(aln[ref_ind,orig_res_ind] + str(ref_res_ind))
            
        aln_res_of_int.append(aln[ref_ind,orig_res_ind] + str(orig_res_ind))
    
    return ref_res_of_int, aln_res_of_int, ref_seq

def write_marked_trees(labeled_trees_fname, tree_node_labels, nodes_to_label, write_file=True):
    #write a file with trees marked
    #Also returns a Dataframe which has node labels next to orthogroup labels

    t_node_label = EvolTree(tree_node_labels, format=8)

    node_names = []
    node_names_paml_in = []
    node_ids_ete = []
    node_ids_paml = []


    for node in t_node_label.traverse():
        node_name = node.name
        node_name_sp = node_name.split('_')
        if len(node_name_sp) == 1: 
            node_id_paml = node_name
            node_name_paml_in = None
        else: 
            node_id_paml = node_name_sp[0]
            node_name_paml_in = node_name_sp[1] + '_' + node_name_sp[2]
        node_names.append(node_name)
        node_ids_ete.append(node.node_id)
        node_ids_paml.append(node_id_paml)
        node_names_paml_in.append(node_name_paml_in)

    node_id_df = pd.DataFrame.from_dict({'name': node_names, 
                                         'name_paml_in': node_names_paml_in, 
                                         'id_paml': node_ids_paml,
                                         'id_ete': node_ids_ete,
                                        }, orient = 'columns')

    id_paml_2_id_ete = dict(zip(node_id_df['id_paml'], node_id_df['id_ete']))
    node_rename_dict = dict(zip(node_id_df['name'], node_id_df['name_paml_in']))

    for node in t_node_label.iter_leaves(): 
        node.name = node_rename_dict[node.name]

    if write_file: 
        with open(labeled_trees_fname,'w') as f_out: 
            f_out.write('  ' + str(len(t_node_label.get_leaf_names())) + '  ' + str(len(nodes_to_label)) + '\n')

            for node in nodes_to_label: 
                t_node_label_marked = t_node_label.copy()
                t_node_label_marked.mark_tree([str(id_paml_2_id_ete[node])],marks = ['#1'])
                f_out.write(t_node_label_marked.write() + '\n')
            
    return node_id_df

def load_tree_with_node_labels(m1_rst_fname):
    #Loads tree with node labels from rst file in ancestral reconstruction of free-branch model DN/DS calculation

    with open(m1_rst_fname,'r') as f_in: 
        line = ''
        while "tree with node labels for Rod Page's TreeView" not in line: 
            line = next(f_in)
        tree_node_labels = next(f_in).strip()

    return tree_node_labels

def load_dn_ds_m1(m1_fname):
    #Loads table of DN DS values for each branch

    with open(m1_fname,'r') as f_in: 
        line = ''
        while "dN & dS for each branch" not in line: 
            line = next(f_in)
        next(f_in)
        headers = next(f_in)
        next(f_in)
        line = next(f_in)
        data = []
        while len(line.split())==9:  #Keep going until there is no more data
            data.append(line.split())
            line = next(f_in)
        data_df = pd.DataFrame(data, columns = headers.split())
        data_df.set_index('branch', inplace=True)
        data_df = data_df.astype('float')

    return data_df

def m1_hits(og_refs, dn_ds_thresh, max_branches_to_test, save_results = True):
    #og_refs: dictionary of gene name and og_refs which have had M1 analysis performed
    #dn_ds_thresh: threshold of dn/ds for a branch to be tested
    #max_branches_to_test: maximum number of branches to test


    marked_tree_branch_info = {}
    hits_tables = {}
    branch_leaves = {}

    for goi, og_ref in og_refs.items(): 
        #etc_og_refs.items(): 
        print(goi)
        #Load tree with node labels
        m1_rst_fname = base_dir + os.sep + os.path.normpath('selection_calculations/m1/' + og_ref + '/rst')
        tree_node_labels = load_tree_with_node_labels(m1_rst_fname)

        #Load table of M1 values
        m1_fname = base_dir + os.sep + os.path.normpath('selection_calculations/m1/' + og_ref + '/m1.out')
        dn_ds_table = load_dn_ds_m1(m1_fname)
        
        #Load map of og names
        og_name_map_fname = base_dir + os.sep + os.path.normpath('msas/structural/tm_align/seq_name_map/' + og_ref + '.tm.tsv')
        og_name_map = pd.read_table(og_name_map_fname)
        
        #Identify Nodes with high DN/DS
        hits_table = dn_ds_table[dn_ds_table['dN/dS']>dn_ds_thresh].sort_values(by='dN/dS', ascending=False)
        print('Total branches above threshold: ' + str(len(hits_table)))
        if len(hits_table)>max_branches_to_test: 
            print('More than ' + str(max_branches_to_test) + ' candidate branches for positive selection, selecting top ' + str(max_branches_to_test) + '. og_ref: ' + og_ref)
            branches_to_test = hits_table.index[0:max_branches_to_test]
        else: 
            branches_to_test = hits_table.index

        hits_table = hits_table.loc[branches_to_test] # truncate hits table
        
        nodes_to_label = [branch.split('..')[1] for branch in hits_table.index]
        
        hits_table['labelled_node'] = nodes_to_label
        
        hits_tables[goi] = hits_table

        #Make file with trees labeling the indicated branches
        labeled_trees_fname = base_dir + os.sep + os.path.normpath('selection_calculations/branch_site/trees/' + og_ref + '.bs_trees')
        marked_tree_branch_info_goi = write_marked_trees(labeled_trees_fname, tree_node_labels, nodes_to_label, write_file = save_results)
        marked_tree_branch_info[goi] = marked_tree_branch_info_goi
        
        branch_leaves_goi = {}
        for jj, selected_node_name in enumerate(nodes_to_label): 
            tree_no = jj + 1
            branch_leaves_goi[tree_no] = get_branch_leaves(tree_node_labels, selected_node_name, og_name_map, marked_tree_branch_info_goi)
        branch_leaves[goi] = branch_leaves_goi


    #Save hits tables in csv files:
    if save_results: 
        for goi, hits_table in hits_tables.items():
            hits_table.to_csv( base_dir + os.sep + os.path.normpath('selection_calculations/branch_site/m1_hits_' + goi + '.csv'))
    
    return hits_tables, branch_leaves, marked_tree_branch_info

def extract_ML_est_BS(bs_fname):
    #output: dictionary of tree_number: (ML values, number beside ML value related to variability of estimate, flag for check convergence issues)

    with open(bs_fname,'r') as f_in: 

        ml_data_out = {}

        for line in f_in:
            if 'TREE' in line: 
                check_conv_flag = False
                tree_no = int(line.split()[2].split(':')[0])
                ml_line = next(f_in)
                if 'check convergence' in ml_line:
                    print('check convergence flag, tree number ' + str(tree_no) + ' in ' + bs_fname)
                    ml_line = next(f_in)
                    check_conv_flag = True
                ml_line_sp = ml_line.split(':')[3].split()
                ml_data_out[tree_no] = (float(ml_line_sp[0]), float(ml_line_sp[1]), check_conv_flag)
    
    return ml_data_out

def parse_branch_site(og_refs, chi2_thresh, beb_val_thresh, p_adj_thresh, hits_tables, branch_leaves): 
    #og_refs = {'POX1': 'OG1122_REF_Scer_AF-P13711-F1-model_v2'}
    #chi2_thresh: Chi2 with 50:50 mixture of point mass 0 has critical values of 5% at 2.71 and 5.41 at 1%.   (page 30 of manual)
    #    For this, a p-value is calculated by putting the test statistic into 1-CDF and dividing by two.  
    #    Normal Chi2 has critical values 3.84 and 5.99 this is what they recommend
    #    p_adj is obtained using stats.false_discovery_control for the p-values obtained for each orthogroup.  
    #
    #beb_val_thresh:  Bayes Empirical value for probability that a residue is positively selected.  Sometimes this value has a high baseline.  
    # 
    #p_adj_thresh: Adjusted p_value threshold
    #
    #Need to initialize branch_leaves and hit_tables with m1_hits
    #
    #Also hits tables are saved in the branch site folder
    # hits_tables = {}
    # for goi, og_ref in og_refs.items():
    #     hits_tables[goi] = pd.read_csv( base_dir + os.sep + os.path.normpath('selection_calculations/branch_site/m1_hits_' + goi + '.csv'))



    bs_analysis_data = {
                'chi2_thresh': chi2_thresh,
                'beb_val_thresh': beb_val_thresh,
                'p_adj_thresh': p_adj_thresh
    }


    for goi, og_ref in og_refs.items(): 
        #goi='QCR2'
        print(goi)
        bs_analysis_data_goi = {}
        #og_ref = etc_og_refs[goi]

        bs_rst_file = base_dir + os.sep + os.path.normpath('selection_calculations/branch_site/' + og_ref + '/A/rst')
        beb_vals = extract_beb_values(bs_rst_file)

        bsa_dir = base_dir + os.sep + os.path.normpath('selection_calculations/branch_site/' + og_ref)

        bsa_fname = bsa_dir + os.sep + os.path.normpath('A/BSA.out')

        ml_data_A = extract_ML_est_BS(bsa_fname)

        bsafix1_fname = bsa_dir + os.sep + os.path.normpath('Afix1/BSAfix1.out')
        ml_data_Afix1 = extract_ML_est_BS(bsafix1_fname)

        p_vals = []
        tree_order = []
        for jj, (tree_no, ml_data_A_tree) in enumerate(ml_data_A.items()):
            tree_order.append(tree_no)
            mlA = ml_data_A_tree[0]
            mlA_check_conv = ml_data_A_tree[2]
            mlAfix1 = ml_data_Afix1[tree_no][0]
            mlAfix1_check_conv = ml_data_Afix1[tree_no][2]
            bs_stat = 2*(mlA-mlAfix1)
            p_val = (1-stats.chi2.cdf(bs_stat,1))/2
            p_vals.append(p_val)
            bs_analysis_data_goi[tree_no] = {'branch_site_stats': 
                                                    {
                                                        'mlA': mlA, 
                                                        'mlA_check_convergence':mlA_check_conv,
                                                        'mlAfix1': mlAfix1,
                                                        'mlAfix1_check_convergence':mlAfix1_check_conv,
                                                        'bs_stat': bs_stat,
                                                        'sig_flag_chi2': bs_stat>chi2_thresh,
                                                        'p_val': p_val 
                                                    }
                                            }

        #get p_adj using p_vals
        p_adj = stats.false_discovery_control(p_vals)

        for jj, tree_no in enumerate(tree_order):
            bs_analysis_data_goi[tree_no]['branch_site_stats']['p_adj'] = p_adj[jj]
            bs_analysis_data_goi[tree_no]['branch_site_stats']['sig_flag_p_adj'] = bool(p_adj[jj]<p_adj_thresh)

        for jj, (branch, row) in enumerate(hits_tables[goi].iterrows()): 
            tree_no = jj+1
            bs_analysis_data_goi[tree_no]['branch_label'] = branch
            bs_analysis_data_goi[tree_no]['labelled_node'] = row['labelled_node']
            bs_analysis_data_goi[tree_no]['branch_site_stats']['m1_dnds'] = row['dN/dS']
            bs_analysis_data_goi[tree_no]['branch_site_stats']['m1_dn'] = row['dN']
            bs_analysis_data_goi[tree_no]['branch_site_stats']['m1_ds'] = row['dS']

        #For each tree add leaf data from bs_analysis_data_goi and extract beb values to identify residues of interest
        for tree_no, branch_leaves_tree in branch_leaves[goi].items():
            bs_analysis_data_goi[tree_no]['leaf_data'] = branch_leaves_tree

            beb_vals_tree = beb_vals[tree_no]
            beb_vals_filt = beb_vals_tree[beb_vals_tree['class_2']>beb_val_thresh]
            residues_of_interest = list(zip(beb_vals_filt.index,beb_vals_filt['ref_AA']))
            pos_sel_res_ref, pos_sel_res_aln, ref_seq = identify_ref_residue(og_ref, residues_of_interest)

            bs_analysis_data_goi[tree_no]['site_data'] = {
                    'pos_sel_res_paml': [ind_res[1] + str(ind_res[0]) for ind_res in residues_of_interest], ####
                    'pos_sel_res_aln': pos_sel_res_aln,
                    'pos_sel_res_ref': pos_sel_res_ref,
                    'beb_score': list(beb_vals_filt['class_2'])
                }



            bs_analysis_data[goi] = bs_analysis_data_goi
    
    return bs_analysis_data

def get_branch_leaves(tree_node_labels, selected_node_name, og_name_map, marked_tree_branch_info_goi):
    #For a given branch site tree and a selected node name, output 
    # leaves_short: a list of proteins (abbreviated name)
    # leaves_long:  a list of proteins (full name)
    # specs_unique: a list of species  
    # paralogs_flag: true if there are paralogs in that branch

    full_tree = Tree(tree_node_labels, format=8)

    for node in full_tree.traverse(): 
        #check if node is a leaf 
        node_name_sp = node.name.split('_')
        if len(node_name_sp)==1:  # node is internal
            node_check = node.name
        else: 
            node_check = node_name_sp[0]

        if selected_node_name == node_check: 
            selected_node = node

    leaves_m1_labels = selected_node.get_leaf_names()
    leaves_short = list(marked_tree_branch_info_goi[marked_tree_branch_info_goi['name'].isin(leaves_m1_labels)]['name_paml_in'])
    leaves_long = list(og_name_map[og_name_map['seq_no'].isin(leaves_short)]['seq_name']) #keeping .pdb label
    specs_non_unique = []
    for seq_name in leaves_long: 
        spec, prot_id = species_from_fasta_id(seq_name.split('.pdb')[0]) # drop .pdb label for this
        specs_non_unique.append(spec)

    specs = list(set(specs_non_unique))

    if len(specs_non_unique)==len(specs):
        paralogs_flag = False
    else: 
        paralogs_flag = True

    branch_leaves_out = {'leaves_short': leaves_short,
                         'leaves_long': leaves_long, 
                         'specs': specs, 
                         'paralogs_flag': paralogs_flag
                        }
    return branch_leaves_out
    
def species_from_fasta_id(fasta_id):
    #Note:  Updated so it no longer takes .pdb. This may break some previous code. 
    #Also updated so that it returns prot_id as well which may break previous code 
    
    #returns the species name given the fasta header (e.g. zygosaccharomyces_rouxii__OG2006__342_4561 or Calb_AF-A0A1D8PDP8-F1-model_v2)
    #Can also handle REF_zygosaccharomyces_rouxii__OG2006__342_4561, although that will only come from clusters with no Scer proteins. 

    fasta_id_sp = fasta_id.split('_')
    if fasta_id_sp[0] in ['REF','Scer','Calb','Spom']:
        if fasta_id_sp[0] == 'REF':
            if fasta_id_sp[1] in ['Scer','Calb','Spom']: 
                spec = model_spec_lookup[fasta_id_sp[1]]
                prot_id = fasta_id_sp[2].split('-')[1]
            else:   
                # ref is not Scer, Calb, Spom
                print('REF seq is not Scer, Calb or Spom')
                spec = '_'.join(fasta_id.split('__')[0].split('_')[1:3])
                prot_id = fasta_id.split('__')[2]
        else:
            spec = model_spec_lookup[fasta_id_sp[0]]
            prot_id = fasta_id_sp[1].split('-')[1]
    else: 
        spec = fasta_id.split('__')[0]
        prot_id = fasta_id.split('__')[2]

    return spec, prot_id

def protein_id_shorten(old_name):
    #Input the name without the .pdb suffix
    
    spec, protein_id = species_from_fasta_id(old_name)
    if spec == 'yHMPu5000034957_hanseniaspora_osmophila_160519': 
        spec = 'hanseniaspora_osmophila'
    elif spec == 'yHMPu5000034604_sporopachydermia_lactativora_160519': 
        spec = 'sporopachydermia_lactativora'
                
    new_name = spec + '__' + protein_id

    return new_name    

def make_model_organism_lookup_from_tm_align_og_ref(save_json=False): 
    #Make Lookup table for model organism genes to og_ref from tm_align clusters
    #Also makes cluster_gene_lists.json which lists all gene ids for all clusters. 
    
    struct_aln_dir = base_dir + os.sep + os.path.normpath('msas/structural/tm_align/fasta_renamed')

    swissprot_id_to_gene_id = load_model_swissprot_id_2_gene_id()

    og_refs =os.listdir(struct_aln_dir)

    struct_align_cluster_gene_lists = {}
    
    SC_orfs_lookup, SC_genename_lookup, SC_features_lookup = read_SGD_features()
    
    for og_ref_fasta in og_refs:
    #og_ref_fasta = og_refs[0]
        og_ref = og_ref_fasta.split('.')[0]

        full_seq_list = []
        scer_seq_list = {'id': [], 'uniprot_id':[], 'genename': [], 'common_name': []}
        #Would be good to have this for C. albicans and S. pombe as well. 

        for record in SeqIO.parse(struct_aln_dir + os.sep + og_ref_fasta, 'fasta'): 
            seq_id = record.id
            full_seq_list.append(seq_id.split('.')[0])
            #Test if sequence is an S.cer sequence.  If so extract the full id, uniprot_id, genename, and common_name

            spec = species_from_fasta_id(seq_id)
            if spec == 'saccharomyces_cerevisiae': 
                scer_seq_list['id'].append(seq_id)
                uniprot_id = seq_id.split('-')[1]
                scer_seq_list['uniprot_id'].append(uniprot_id)
                genename= swissprot_id_to_gene_id['Scer'][uniprot_id]
                scer_seq_list['genename'].append(genename)
                scer_seq_list['common_name'].append(SC_genename_lookup[genename])


        struct_align_cluster_gene_lists[og_ref] = {'all_sequences': full_seq_list, 'scer_sequences': scer_seq_list} 

    struct_align_cluster_gene_lists_fname = base_dir + os.sep + os.path.normpath('msas/structural/tm_align/cluster_gene_lists.json')
    if save_json: 
        with open(struct_align_cluster_gene_lists_fname, 'w') as f:
            json.dump(struct_align_cluster_gene_lists, f, sort_keys=True, indent=4 )
    
    return struct_align_cluster_gene_lists

def surface_core_analysis_untrimmed_to_trimmed(aln_trim_fname,surface_core_aligned_mat_untrimmed):
    zerocolumns = False
    aln_trimmed = AlignIO.read(open(aln_trim_fname),'fasta')
    
    #Load clipkit log
    aln_trim_log_fname = aln_trim_fname + '.log'
    aln_trim = AlignIO.read(open(aln_trim_fname),'fasta')
    
    record_order = []
    for record in aln_trim:
        record_order.append(record.id)
        
    aln_trim_log = pd.read_table(aln_trim_log_fname, sep=' ', header=None, index_col=0)
    aln_trim_log.rename(columns={1:'trim',2:'note',3:'value'}, inplace=True)

    #makes a zero-based list of indices of the original full alignment index for extracting the trimmed index
    #trim_ind_2_orig_ind = {}
    aln_trim_log_kept = aln_trim_log[aln_trim_log['trim']=='keep']
    aln_trim_log_kept_list = []
    for trim_ind, orig_ind in enumerate(aln_trim_log_kept.index): 
        aln_trim_log_kept_list.append(orig_ind-1)
        #trim_ind_2_orig_ind[trim_ind+1] = orig_ind

    surface_core_aligned_mat = surface_core_aligned_mat_untrimmed[:,aln_trim_log_kept_list]
    
    core_column_pct = [] #Surface column pct will be 1-core column pct
    column_occupancy = []
    nspecs = surface_core_aligned_mat.shape[0]
    for jj in range(surface_core_aligned_mat.shape[1]): 
        col = surface_core_aligned_mat[:,jj]
        col_count = Counter(col)
        if (col_count['C']+col_count['S'])==0: 
            zerocolumns = True  
            core_column_pct.append(np.nan)
        else: 
            core_column_pct.append(col_count['C']/(col_count['C']+col_count['S']))
        column_occupancy.append((nspecs-col_count['-'])/nspecs)
    
    if zerocolumns: 
        print('at least one trimmed alignment column has no surface or core residues - may be empty')
        
    return core_column_pct, column_occupancy, surface_core_aligned_mat, record_order

def surface_core_analysis_alignment(aln, og_surface_core_data):
    #Given an untrimmed alignment, extracts the percentage of residues that are core, and the percentage of columns of the alignment that are occupied, and the full matrix with "I", "O", or "-" in each column.  
    #Also provides the order of the proteins from the alignment which may not be the same as in the core/surface data
    zerocolumns = False
    
    surface_core_aligned = []
    record_order = []
    for record in aln:
        record_order.append(record.id)
        seq_out, msa2input, input2msa, pair_mapping = seq_squeeze(record.seq)

        surface_core_unaligned_prot = og_surface_core_data['protein_data'][record.id.split('.')[0]]['full_protein_assignment']

        surface_core_aligned_prot = []

        for jj, res in enumerate(record.seq): 
            if res=='-': 
                surface_core_aligned_prot.append('-')
            else: 
                surface_core_aligned_prot.append(surface_core_unaligned_prot[msa2input[jj]])

        surface_core_aligned.append(surface_core_aligned_prot)

    surface_core_aligned_mat = np.array(surface_core_aligned)

    core_column_pct = [] #Surface column pct will be 1-core column pct
    column_occupancy = []
    nspecs = surface_core_aligned_mat.shape[0]
    for jj in range(surface_core_aligned_mat.shape[1]): 
        col = surface_core_aligned_mat[:,jj]
        col_count = Counter(col)
        if (col_count['C']+col_count['S'])==0: 
            zerocolumns=True
            core_column_pct.append(np.nan)
        else: 
            core_column_pct.append(col_count['C']/(col_count['C']+col_count['S']))
        column_occupancy.append((nspecs-col_count['-'])/nspecs)
    
    if zerocolumns: 
        print('alignment column has no surface or core residues - may be empty') 
    
    return core_column_pct, column_occupancy, surface_core_aligned_mat, record_order

def extract_surface_core(og_summary, sasa_cut, dict_sasa_max):
    #For a given protein lists how many are Core, Surface, Missing From Alignment, Not assigned 
    #Would be good to check Oliver's alignment vs my alignment to make sure it is correct.  
    #This should be something saved somewhere as a reference
    
    log = ''

    og_surface_core_data = {}
    og_surface_core_data['protein_data']={}

    for prot_ind, prot_name in enumerate(og_summary['Molecules'].keys()): 
        #prot_name = 'yHMPu5000034957_hanseniaspora_osmophila_160519__OG1316__247_2096'
        #prot_ind = prot_ind_lookup[prot_name]
        mol_sasa = og_summary['Molecules'][prot_name]['SASA']
        mol_res = og_summary['Molecules'][prot_name]['Residues']

        core_res_all = (mol_sasa/np.vectorize(dict_sasa_max.get)(mol_res))<sasa_cut
        surface_res_all = (mol_sasa/np.vectorize(dict_sasa_max.get)(mol_res))>sasa_cut

        #Initialize surface_core, a vector that indicates whether a residue in the sequence is, core, surface, or NA if not calculated 
        surface_core = ['NA']*len(mol_res)

        #surface_residues: 
        for jj, surface_bool in enumerate(surface_res_all): 
            if surface_bool: 
                #seq_ind = mapping_ref_based_aln[jj]
                surface_core[jj]= 'S'

        #core_residues: 
        for jj, core_bool in enumerate(core_res_all): 
            if core_bool: 
                #seq_ind = mapping_ref_based_aln[jj]
                surface_core[jj]= 'C'

        surface_core_count = Counter(surface_core)
        
        missing_core_surface = {'C','S'} - set(surface_core_count.keys())
        
        #check to see if either all Core or all Surface residues are missing:  
        if len(missing_core_surface)>0: 
            message = prot_name + ' missing ' + ' and '.join([feature_type for feature_type in list(missing_core_surface)])
            print(message)
            log = log + message + '\n'

            for c_s in list(missing_core_surface):
                surface_core_count[c_s] = 0

        #Check that all residues for the protein in the ref_based alignment are assigned to either core or surface. 
        assert len(mol_res) == (surface_core_count['S'] + surface_core_count['C']), 'Residue aligned to reference structure not assigned as core or surface'

        
        surface_core_pct = {}
        for key, count in surface_core_count.items():
            surface_core_pct[key] = count/len(mol_res)

             #return surface_core_pct, surface_core_count, surface_core

        og_surface_core_data['protein_data'][prot_name] = {'full_protein_assignment' : surface_core,
                                          'full_protein_count': {'core': surface_core_count['C'], 'surface': surface_core_count['S']}, 
                                          'full_protein_pct' :  {'core': surface_core_pct['C'], 'surface': surface_core_pct['S']}
                                         }
    
    cores = []
    surfaces = []

    for key, data in og_surface_core_data['protein_data'].items(): 
        #print(key)
        cores.append(data['full_protein_pct']['core'])
        surfaces.append(data['full_protein_pct']['surface'])

    og_surface_core_data['core_mean'] = np.mean(cores)
    og_surface_core_data['core_std'] = np.std(cores)
        
    return og_surface_core_data, log

def load_hyphy_busted_data(hyphy_base_dir):
    # Extract Hyphy busted data from a directory base containing outputs

    og_list = []
    og_ref_list = []
    hyphy_data_out = {}
    partitions = ['0','1','2']
    partition_label = {'0':'low', '1':'med', '2': 'high'}

    for og_ref in next(os.walk(hyphy_base_dir))[1]:    
        #next(os.walk(m0_base))[0] gets just directories - https://stackoverflow.com/questions/141291/how-to-list-only-top-level-directories-in-python
        #og_ref = 'OG1299_REF_Scer_AF-P00549-F1-model_v2'
        og = og_ref.split('_')[0]
        og_ref_list.append(og_ref)
        og_list.append(og)
        #print(og_ref)
        hyphy_dir = hyphy_base_dir + og_ref + os.sep

        hyphy_json = hyphy_dir + 'hyphy_busted.json'
        
        try: 
            with open(hyphy_json, 'r') as f:
                hyphy_data = json.load(f) 

            average_omega = 0
            average_rate = 0

            unconstrained_model_data = hyphy_data['fits'][ 'Unconstrained model']['Rate Distributions']

            output_data = {}
            for partition in partitions: 
                output_data['omega_' + partition_label[partition]] = unconstrained_model_data['Test'][partition]['omega']
                output_data['omega_' + partition_label[partition] + '_proportion'] = unconstrained_model_data['Test'][partition]['proportion']
                output_data['rate_' + partition_label[partition]] = unconstrained_model_data['Synonymous site-to-site rates'][partition]['rate']
                output_data['rate_' + partition_label[partition] + '_proportion'] = unconstrained_model_data['Synonymous site-to-site rates'][partition]['proportion']
                output_data['pvalue'] = hyphy_data['test results']['p-value']
                average_omega = average_omega + output_data['omega_' + partition_label[partition]]*output_data['omega_' + partition_label[partition] + '_proportion']
                average_rate = average_rate + output_data['rate_' + partition_label[partition]]*output_data['rate_' + partition_label[partition] + '_proportion']

            output_data['omega_avg']=average_omega
            output_data['rate_avg']= average_rate

            hyphy_data_out[og_ref] = output_data
        
        except FileNotFoundError:
            print(og_ref + ' has no hyphy BUSTED output')
            

    hyphy_data_df = pd.DataFrame.from_dict(hyphy_data_out, orient='index')

    return hyphy_data_df