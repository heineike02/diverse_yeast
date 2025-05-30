#Remake alignments after filtering sparse sequences to test M0

#import sys
import os
import subprocess
from Bio import SeqIO
from ete3 import Tree
import pandas as pd
import shutil
import numpy as np

base_dir = os.path.normpath('/home/heineike_wsl2/alphafold')
#aln_dir = base_dir + os.sep + os.path.normpath('msas/structural/tm_align/fasta_renamed') 
#align_files = os.listdir(aln_dir)

selected_alignments = ['OG1273_REF_Scer_AF-P40012-F1-model_v2','OG3104_REF_Scer_AF-P53954-F1-model_v2','OG1230_REF_Scer_AF-P29465-F1-model_v2', 'OG1364_REF_Scer_AF-P27796-F1-model_v2']
#[fname.split('.')[0] for fname in align_files]

seq_filt_thresh = 0.85 #Threshold percentage of residues that must be present

#trim_msa_thresh = 0.25  # Threshold to remove clusters that have poor alignments.  If the strict trimming MSA length is less than .25 * median sequence length, the cluster is removed. 
#ogs_filtered = []

#os.path.normpath('/home/heineike_wsl2/Crick_LMS/projects/diverse_yeasts/alphafold')
# output_dir = base_dir + os.sep + os.path.normpath('selection_calculations/20220526_sel_calc')

#selected_ogs = ['OG2645']#['OG4150', 'OG2603', 'OG3677', 'OG2845']

#selected_og_refs = ['OG2645_REF_Scer_AF-P05375-F1-model_v2']#['OG4150_REF_Scer_AF-P07256-F1-model_v2', 'OG2603_REF_Scer_AF-P50076-F1-model_v2', 'OG2845_REF_Scer_AF-P43577-F1-model_v2', 'OG3677_REF_Scer_AF-P47125-F1-model_v2', 'OG1299_REF_Scer_AF-P00549-F1-model_v2']


## Get this to handle non REF ones or filter them out
#OG1111_alloascoidea_hylecoeti__OG1111__0_3867
#og_ref = 'OG4150_REF_Scer_AF-P07256-F1-model_v2'

for alignment in selected_alignments: 
    print(alignment)
    #     if alignment.split('_')[1]=='REF':
    #         og,ref = alignment.split('_REF_')
    #     else: 
    #     og = alignment.split('_')[0]

    #     og_pep_msa_fname = base_dir + os.sep + os.path.normpath('msas/structural/tm_align/fasta_renamed/' + alignment + '.tm.fasta')

    #     #strict trimming for codon alignments    
    #     print('Strict trimming')
    #     clipkit_cmd = ['clipkit', og_pep_msa_fname, '-m','gappy' ,'-g','0.1', '-l']
    #     subprocess.run(clipkit_cmd)

    #     #Move clipkit log and output
    #     suffixes = ['.clipkit','.clipkit.log' ]

    #     for suffix in suffixes: 
    #         #Move clipkit files to new folder
    #         shutil.move(og_pep_msa_fname + suffix, 
    #                     base_dir + os.sep + os.path.normpath('msas/structural/tm_align/trim_strict/' + alignment + '.tm.fasta' + suffix)
    #                    )

    #     ##Open alignment and check the length.  If the length is less than 25% of the length of the ref structure then skip and throw a flag.  
    #     msa_pep_trimmed = base_dir + os.sep + os.path.normpath('msas/structural/tm_align/trim_strict/' + alignment + '.tm.fasta.clipkit')

    #     og_aln_fasta = SeqIO.parse(og_pep_msa_fname, 'fasta')
    #     msa_pep_trimmed_fasta = SeqIO.parse(msa_pep_trimmed,'fasta')

    #     lengths = []
    #     for record in og_aln_fasta: 
    #         lengths.append(int(record.description.split('\t')[1].split('=')[1]))
    #     med_length = np.median(lengths)

    #     first_rec = next(msa_pep_trimmed_fasta)
    #     trimmed_msa_len = len(first_rec.seq)

    #     if (trimmed_msa_len/med_length)>trim_msa_thresh:    
    #         #Thread original alignment               
    #         ##Good place for a test that these alignments have same AA sequence as the protein sequence and CDS.  

    #         print('Thread orig alignment')
    #         og_cds_fname = base_dir + os.sep +  os.path.normpath('selected_proteins/og_sequences/cds_tm/' + alignment + '.cds.fasta')
    #         og_cds_msa_fname = base_dir + os.sep +  os.path.normpath('msas/structural/tm_align/cds_aln/' + alignment +  '.tm.cds.aln.fasta')

    #         phykit_cmd = ['phykit', 'thread_dna',
    #                       '-p', og_pep_msa_fname,
    #                       '-n', og_cds_fname, 
    #                      ]

    #         with open(og_cds_msa_fname,'w') as f_cds:
    #             output = subprocess.run(phykit_cmd, stdout=f_cds)               

    #         #Make Trimmed Alignment
    #         print('Make Trimmed Alignment')

    #         msa_cds_trimmed = base_dir + os.sep +  os.path.normpath('msas/structural/tm_align/cds_trim_strict/' + alignment +  '.tm.fasta.clipkit.cds')

    #         phykit_cmd = ['phykit', 'thread_dna',
    #                       '-p', msa_pep_trimmed,
    #                       '-n', og_cds_msa_fname, 
    #                       '-c', msa_pep_trimmed + '.log'
    #                      ]

    #         with open(msa_cds_trimmed,'w') as f_cds_trimmed:
    #             output = subprocess.run(phykit_cmd, stdout=f_cds_trimmed)


    #         ##Rename cds alignment 
    #         print('Rename')
    #         msa_cds_trimmed_renamed = msa_cds_trimmed+'.renamed'

    #         #shorten the name of the fasta and build name_replace dictionary

    #         #Make the sequence name map
    #         og_pep_msa = SeqIO.parse(og_pep_msa_fname,'fasta')
    #         seq_name_map_fname = base_dir + os.sep +  os.path.normpath('msas/structural/tm_align/seq_name_map/' + alignment + '.tm.tsv')
    #         with open(seq_name_map_fname,'w') as fout_seq_name: 
    #             fout_seq_name.write('seq_name\tseq_no\n')
    #             for ind, record in enumerate(og_pep_msa):
    #                 seq_id = record.id.split('.')[0]
    #                 fout_seq_name.write(record.id + '\t' + og + '_'+ str(ind)+'\n')

    #         seq_name_map_df = pd.read_csv(seq_name_map_fname, sep='\t')
    #         seq_name_map = dict(zip(seq_name_map_df['seq_name'],seq_name_map_df['seq_no']))

    #         msa_cds_trimmed_seqs = SeqIO.parse(msa_cds_trimmed, 'fasta')
    #         with open(msa_cds_trimmed_renamed,'w') as f_out: 
    #             for record in msa_cds_trimmed_seqs:
    #                 new_id = seq_name_map[record.id] 
    #                 f_out.write('>' + str(new_id) + '\n')
    #                 f_out.write(str(record.seq) + '\n')
    

    #remove sequences that have many missing residues
    #Assume trimmed sequence already created 
    msa_cds_trimmed = base_dir + os.sep +  os.path.normpath('msas/structural/tm_align/cds_trim_strict/' + alignment +  '.tm.fasta.clipkit.cds')
    msa_cds_trimmed_renamed = msa_cds_trimmed+'.renamed' 
    msa_cds_trimmed_seqs = SeqIO.parse(msa_cds_trimmed_renamed, 'fasta')

    seq_ids = []
    aln_len_sq = [] 

    for record in msa_cds_trimmed_seqs: 
        seq_ids.append(record.id)
        seq_squeezed = len([res for res in record.seq if res!='-'])
        aln_len_sq.append(seq_squeezed)
    aln_len = len(record.seq)

    aln_len_df = pd.DataFrame.from_dict({'seq_id': seq_ids,'seq_len_squeeze':aln_len_sq})
    aln_len_df['seq_len_pct'] = aln_len_df['seq_len_squeeze']/aln_len

    aln_len_df_filt = aln_len_df[aln_len_df['seq_len_pct']>seq_filt_thresh]

    msa_cds_trimmed_seqs = SeqIO.parse(msa_cds_trimmed_renamed, 'fasta')
    msa_cds_trimmed_filtered = base_dir + os.sep +  os.path.normpath('selection_calculations/m0_yn00_comparison/m0/alignments/' + alignment +  '.m0_filtered.fasta')

    with open(msa_cds_trimmed_filtered, 'w') as f_out:
        for record in msa_cds_trimmed_seqs: 
            seq_id = record.id
            if seq_id in set(aln_len_df_filt['seq_id']):
                print(seq_id)
                f_out.write('>' + seq_id + '\n')
                f_out.write(str(record.seq) + '\n')
    
       
    print('Convert to Phylip')
    #converts fasta to PHYLIP format 
    msa_cds_trimmed_phy = base_dir + os.sep +  os.path.normpath('selection_calculations/m0_yn00_comparison/m0/alignments/' + alignment +  '.m0_filtered.phy')
    biokit_cmd = ['biokit','file_format_converter',
                  '-i',msa_cds_trimmed_filtered, '-iff', 'fasta',
                  '-o', msa_cds_trimmed_phy,
                  '-off', 'phylip']

    subprocess.run(biokit_cmd)

    #adds I in phylip to indicate interleaved status
    msa_cds_trimmed_phy_codeML = base_dir + os.sep + os.path.normpath('selection_calculations/m0_yn00_comparison/m0/alignments/' + alignment +  '.m0_filtered.codeML.phy')

    with open(msa_cds_trimmed_phy,'r') as f_in:
        with open(msa_cds_trimmed_phy_codeML,'w') as f_out: 
            line = next(f_in)
            line_out = line.strip('\n') + ' I\n'
            f_out.write(line_out)
            for line in f_in:
                f_out.write(line)
    
    
    
    ## Also filter relevant tree
    tree_orig = base_dir + os.sep + os.path.normpath('msas/structural/tm_align/trees/' + alignment + '.tm.fasta.clipkit.treefile')
    tree_renamed = tree_orig+'.renamed'
    tree_filtered = base_dir + os.sep +  os.path.normpath('selection_calculations/m0_yn00_comparison/m0/alignments/' + alignment + '.filtered')
    
    full_tree = Tree(tree_renamed, format=1) #Using flexible with internal node names as support is not in the right format to load with format = 0
    full_tree.prune(list(aln_len_df_filt['seq_id']))
    full_tree.write(format=5, outfile = tree_filtered)
    
    #     else: 
    #         print(alignment + ' skipped: Trimmed alignment shorter than ' + str(trim_msa_thresh) + 'times average alginment length')
    #         ogs_filtered.append(alignment)


# filtered_og_log_fname = base_dir + os.sep +  os.path.normpath('msas/structural/tm_align/cds_trim_strict/filtered_og.log')

# with open(filtered_og_log_fname, 'w') as f_out: 
#     f_out.write('Orthogroups filtered out when running 20221206_struct_align_dnds_msas.py because strict trimming of alignment was below trim_msa_thresh=' + str(trim_msa_thresh) + ' * average sequence length threshold\n')
#     for og in ogs_filtered: 
#         f_out.write(og + '\n')



