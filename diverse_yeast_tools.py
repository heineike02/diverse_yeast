import os
import pandas as pd
import json
import pickle
from Bio import SeqIO
from Bio import AlignIO   #, Align

base_dir = os.path.normpath('G:/My Drive/Crick_LMS/projects/diverse_yeasts/alphafold')
divyeast_dir = os.path.normpath('C:/Users/heineib/Documents/GitHub/diverse_yeast')
y1000plus_dir = os.path.normpath('C:/Users/heineib/Documents/GitHub/y1000plus_tools/data') + os.sep
genomes_dir = os.path.normpath('G:/My Drive/Crick_LMS/external_data/genomes')


#Function to parse sequence id from first line of fasta file in either shen or uniprot proteomes. 

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

    spec_lookup_fname = y1000plus_dir + os.path.normpath('y1000plus_tools_data/y1000plus/id_lookups/' + spec_orig_genome + '.csv')
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
                   'Spom': genomes_dir +os.sep +  os.path.normpath('schizosaccharomyces_pombe\peptide.fa')
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
    scer_lookup_fname = y1000plus_dir + os.path.normpath('y1000plus_tools_data/y1000plus/id_lookups/saccharomyces_cerevisiae.csv')
    scer_lookup = pd.read_csv(scer_lookup_fname, index_col=0)
    gene_id_2_y1000_id['Scer'] = dict(zip(scer_lookup.index,scer_lookup['y1000_id']))

    #Load C.alb lookup table
    calb_lookup_fname = y1000plus_dir + os.path.normpath('y1000plus_tools_data/y1000plus/id_lookups/candida_albicans.csv')
    calb_lookup = pd.read_csv(calb_lookup_fname, index_col=0)
    gene_id_2_y1000_id['Calb'] = dict(zip(calb_lookup.index,calb_lookup['y1000_id']))

    #Spom is difficult because the y1000 id map goes from gene names rather than systematic ids and many of the gene names are synonyms rather than the official names
    #Load S.pom lookup table: 
    spom_lookup_fname = y1000plus_dir + os.path.normpath('y1000plus_tools_data/y1000plus/id_lookups/schizosaccharomyces_pombe.csv')
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
    
    
    #     ###Should update this one
    #     calb_swissprot_id_2_gene_id = pickle.load(open(base_dir + os.sep + os.path.normpath('msas/structural/Mapping_calb.pkl'),"rb"))
    #     swissprot_id_2_gene_id['Calb'] = dict(zip(calb_swissprot_id_2_gene_id.values(),calb_swissprot_id_2_gene_id.keys()))

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
    #Lookup from swissprot_id to gene_id
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

def seq_squeeze(seq_in):
    seq_out = ''
    for res in seq_in:
        if res!= '-':
            seq_out = seq_out + res
    
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