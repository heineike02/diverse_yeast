#Make Strictly trimmed CDS for DN/DS calculation for bundung site only based on subset from tm_align msas
#
#
#Input: 
# Subset alignment and CDS alignment from msas/structural/tm_align/fasta_renamed
# ?Alignments from msas/structural/tm_align/fasta_renamed/feature_subsets/binding_site/
#                       cds_aln
#                       fasta_renamed
#
#Output: 
#  Strict trimming: msas/structural/tm_align/trim_strict/
#  threaded original alignment: msas/structural/tm_align/cds_aln/
#  Trimmed CDS alignment: msas/structural/tm_align/cds_trim_strict/
#  Renaming Map:  msas/structural/tm_align/seq_name_map/
#
#  When running the first time need to make the trim_strict, cds_trim_strict and seq_name_map folders. 
#
#  The final output are phylip files with short names strictly trimmed and formatted for use with Code ML.
#
#  msas/structural/tm_align/cds_trim_strict/<alignment>.tm.fasta.clipkit.cds.renamed.codeML.phy'


#import sys
import os
import subprocess
from Bio import SeqIO
import pandas as pd
import shutil
import numpy as np
#from diverse_yeast_tools import seq_squeeze

def seq_squeeze(seq_in): 
    seq_out = ''
    for res in seq_in:
        if res!= '-':
            seq_out = seq_out + res

    return seq_out
            
base_dir = os.path.normpath('/home/heineikeb/alphafold')
aln_dir = base_dir + os.sep + os.path.normpath('msas/structural/tm_align') 

feature = "uniprot_binding_site"
feature_subset_dir = aln_dir + os.sep + os.path.normpath('feature_subsets/' + feature) 


align_files = os.listdir(feature_subset_dir + os.sep + os.path.normpath('fasta_renamed'))
selected_alignments = [fname.split('.')[0] for fname in align_files]

trim_msa_thresh = 0.25  # Threshold to remove clusters that have poor alignments.  If the strict trimming MSA length is less than .25 * median sequence length, the cluster is removed. 

# As of 10 Jan 2023 that resulted in three OGs being removed. 
# Orthogroups filtered out when running 20221206_struct_align_dnds_msas.py because strict trimming of alignment was below trim_msa_thresh=0.25 * average sequence length threshold
# OG2147_REF_Scer_AF-P39692-F1-model_v2
# OG1306_REF_Scer_AF-P38298-F1-model_v2
# OG1746_REF_Scer_AF-P32642-F1-model_v2


ogs_filtered = []

#os.path.normpath('/home/heineike_wsl2/Crick_LMS/projects/diverse_yeasts/alphafold')
# output_dir = base_dir + os.sep + os.path.normpath('selection_calculations/20220526_sel_calc')

#selected_ogs = ['OG2645']#['OG4150', 'OG2603', 'OG3677', 'OG2845']

#selected_og_refs = ['OG2645_REF_Scer_AF-P05375-F1-model_v2']#['OG4150_REF_Scer_AF-P07256-F1-model_v2', 'OG2603_REF_Scer_AF-P50076-F1-model_v2', 'OG2845_REF_Scer_AF-P43577-F1-model_v2', 'OG3677_REF_Scer_AF-P47125-F1-model_v2', 'OG1299_REF_Scer_AF-P00549-F1-model_v2']


## Get this to handle non REF ones or filter them out
#OG1111_alloascoidea_hylecoeti__OG1111__0_3867
#og_ref = 'OG4150_REF_Scer_AF-P07256-F1-model_v2'

for alignment in selected_alignments: 
    #     print(alignment)
    #     if alignment.split('_')[1]=='REF':
    #         og,ref = alignment.split('_REF_')
    #     else: 
    og = alignment.split('_')[0]
    
    og_pep_msa_fname = feature_subset_dir + os.sep + os.path.normpath('fasta_renamed/' + alignment + '.tm.fasta')
    
    #strict trimming for codon alignments    
    print('Strict trimming')
    clipkit_cmd = ['clipkit', og_pep_msa_fname, '-m','gappy' ,'-g','0.1', '-l']
    subprocess.run(clipkit_cmd)

    #Move clipkit log and output
    suffixes = ['.clipkit','.clipkit.log' ]

    for suffix in suffixes: 
        #Move clipkit files to new folder
        shutil.move(og_pep_msa_fname + suffix, 
                    feature_subset_dir + os.sep + os.path.normpath('trim_strict/' + alignment + '.tm.fasta' + suffix)
                   )

    ##Open alignment and check the length.  If the length is less than 25% of the length of the ref structure then skip and throw a flag.  
    msa_pep_trimmed = feature_subset_dir + os.sep + os.path.normpath('trim_strict/' + alignment + '.tm.fasta.clipkit')

    og_aln_fasta = SeqIO.parse(og_pep_msa_fname, 'fasta')
    msa_pep_trimmed_fasta = SeqIO.parse(msa_pep_trimmed,'fasta')

    lengths = []
    for record in og_aln_fasta: 
        lengths.append(len(seq_squeeze(str(record.seq))))
        #lengths.append(int(record.description.split('\t')[1].split('=')[1]))
    med_length = np.median(lengths)

    first_rec = next(msa_pep_trimmed_fasta)
    trimmed_msa_len = len(first_rec.seq)

    if (trimmed_msa_len/med_length)>trim_msa_thresh:
# Not needed as we use the threaded alignment to get a CDS subset
#         #Thread original alignment               
#         ##Good place for a test that these alignments have same AA sequence as the protein sequence and CDS.  

#         print('Thread orig alignment')
#         og_cds_fname = base_dir + os.sep +  os.path.normpath('selected_proteins/og_sequences/cds_tm/' + alignment + '.cds.fasta')
#         og_cds_msa_fname = aln_dir + os.sep + os.path.normpath('cds_aln/' + alignment +  '.tm.cds.aln.fasta')

#         phykit_cmd = ['phykit', 'thread_dna',
#                       '-p', og_pep_msa_fname,
#                       '-n', og_cds_fname, 
#                      ]

#         with open(og_cds_msa_fname,'w') as f_cds:
#             output = subprocess.run(phykit_cmd, stdout=f_cds)               

#         #Make Trimmed Alignment
        print('Make Trimmed CDS Alignment')
    
        og_cds_msa_fname = feature_subset_dir + os.sep + os.path.normpath('cds_aln/' + alignment +  '.tm.cds.aln.fasta')
    
        msa_cds_trimmed = feature_subset_dir + os.sep + os.path.normpath('cds_trim_strict/' + alignment +  '.tm.fasta.clipkit.cds')

        phykit_cmd = ['phykit', 'thread_dna',
                      '-p', msa_pep_trimmed,
                      '-n', og_cds_msa_fname, 
                      '-c', msa_pep_trimmed + '.log'
                     ]
        
        print(' '.join(phykit_cmd))
        with open(msa_cds_trimmed,'w') as f_cds_trimmed:
            output = subprocess.run(phykit_cmd, stdout=f_cds_trimmed)


        ##Rename cds alignment 
        print('Rename cds alignment for CodeML')
        msa_cds_trimmed_renamed = msa_cds_trimmed+'.renamed'

        #shorten the name of the fasta and build name_replace dictionary

        #Load the sequence name map
        seq_name_map_fname = aln_dir + os.sep + os.path.normpath('seq_name_map/' + alignment + '.tm.tsv')
        seq_name_map_df = pd.read_csv(seq_name_map_fname, sep='\t')
        seq_name_map = dict(zip(seq_name_map_df['seq_name'],seq_name_map_df['seq_no']))

        msa_cds_trimmed_seqs = SeqIO.parse(msa_cds_trimmed, 'fasta')
        with open(msa_cds_trimmed_renamed,'w') as f_out: 
            for record in msa_cds_trimmed_seqs:
                new_id = seq_name_map[record.id] 
                f_out.write('>' + str(new_id) + '\n')
                f_out.write(str(record.seq) + '\n')

        print('Convert to Phylip')
        #converts fasta to PHYLIP format 
        msa_cds_trimmed_phy = msa_cds_trimmed + '.renamed.phy'
        biokit_cmd = ['biokit','file_format_converter',
                      '-i',msa_cds_trimmed_renamed, '-iff', 'fasta',
                      '-o', msa_cds_trimmed_phy,
                      '-off', 'phylip']

        subprocess.run(biokit_cmd)

        #adds I in phylip to indicate interleaved status
        msa_cds_trimmed_phy_codeML = msa_cds_trimmed + '.renamed.codeML.phy'

        with open(msa_cds_trimmed_phy,'r') as f_in:
            with open(msa_cds_trimmed_phy_codeML,'w') as f_out: 
                line = next(f_in)
                line_out = line.strip('\n') + ' I\n'
                f_out.write(line_out)
                for line in f_in:
                    f_out.write(line)
    
    else: 
        print(alignment + ' skipped: Trimmed alignment shorter than ' + str(trim_msa_thresh) + 'times average alginment length')
        ogs_filtered.append(alignment)
        

filtered_og_log_fname = feature_subset_dir + os.sep + os.path.normpath('cds_trim_strict/filtered_og.log')

with open(filtered_og_log_fname, 'w') as f_out: 
    f_out.write('Orthogroups filtered out when running 20221206_struct_align_dnds_msas.py because strict trimming of alignment was below trim_msa_thresh=' + str(trim_msa_thresh) + ' * average sequence length threshold\n')
    for og in ogs_filtered: 
        f_out.write(og + '\n')



