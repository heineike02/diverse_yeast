#Using the same sequences as the structural alignment, make sequence based alignments using mafft and clustalo. 
#Then: 
#   1. Trim the alignments with standard parameters and assemble trees
#   2. Make a strict trimming and prepare files for DN/DS calculations


import os
import subprocess
from Bio import SeqIO
#import pandas as pd
#import shutil
#import numpy as np

home_dir = os.path.normpath('/home/heineike_wsl2')

base_dir = home_dir + os.sep + 'alphafold'
proteome_dir = base_dir + os.sep + os.path.normpath('selected_proteins/og_sequences/proteome_tm')

proteome_files = os.listdir(proteome_dir)
selected_proteomes = [fname.split('.')[0] for fname in proteome_files]


aln_dir = base_dir + os.sep + os.path.normpath('msas/sequence') 


#for proteome in selected_proteomes: 

proteome = 'OG1254_REF_Scer_AF-P40395-F1-model_v2'  #selected_proteomes[0]

#Run Mafft on proteome with auto parameter.  

#Output the screen as a log so that I can extract the strategy as .aln.fasta.log
# mafft_command = ['mafft', '--auto', 
#                  proteome_dir + os.sep + proteome + '.pep.fasta']

# mafft_output_fname =  aln_dir + os.sep + os.path.normpath('mafft/fasta/'+ proteome + '.aln.fasta')
# mafft_log_fname = aln_dir + os.sep + os.path.normpath('mafft/fasta/'+ proteome + '.aln.fasta.log')

# with open(mafft_output_fname, 'w') as mafft_output: 
#     with open(mafft_log_fname, 'w') as mafft_log: 
#         subprocess.run(mafft_command, stdout=mafft_output, stderr=mafft_log)
    

    
#Run Clustalo on proteome with default settings
clustalo_output_fname =  aln_dir + os.sep + os.path.normpath('clustalo/fasta/'+ proteome + '.aln.fasta')
clustalo_log_fname = aln_dir + os.sep + os.path.normpath('clustalo/fasta/'+ proteome + '.aln.fasta.log')

clustalo_command = [home_dir + os.sep + os.path.normpath('clustal/clustalo-1.2.4-Ubuntu-x86_64'), 
                    '-i', proteome_dir + os.sep + proteome + '.pep.fasta',
                    '-o', clustalo_output_fname,
                    '-v','-v',
                    '-l', clustalo_log_fname
                   ]
subprocess.run(clustalo_command)

print([' '.join(clustalo_command)])


#Therefore I am making a separate bash script and running that 

# mafft_cmd_file_fname = '/home/heineike_wsl2/alphafold/msas/sequence/mafft/mafft_cmd.sh'

# with open(mafft_cmd_file_fname, 'w') as mafft_cmd_file:
#     mafft_cmd_file.write('#!/bin/bash\n')
#     mafft_cmd_file.write(' '.join(mafft_command))

# subprocess.run(mafft_cmd_file_fname)
#'--genafpair', '--maxiterate', '1000', 





# ls_cmd = ['ls', proteome_dir + os.sep]
# subprocess.run(ls_cmd)


#/home/heineike_wsl2/alphafold/selected_proteins/og_sequences/proteome_tm/OG1004_REF_Scer_AF-P15938-F1-model_v2.pep.fasta

# mafft --auto /home/heineike_wsl2/alphafold/selected_proteins/og_sequences/proteome_tm/OG1004_REF_Scer_AF-P15938-F1-model_v2.pep.fasta > /home/heineike_wsl2/alphafold/msas/sequence/mafft/fasta/OG1004_REF_Scer_AF-P15938-F1-model_v2.aln.fasta

# mafft --auto /home/heineike_wsl2/alphafold/selected_proteins/og_sequences/proteome_tm/OG1004_REF_Scer_AF-P15938-F1-model_v2.pep.fasta > /home/heineike_wsl2/alphafold/msas/sequence/mafft/fasta/OG1004_REF_Scer_AF-P15938-F1-model_v2.aln.fasta

#mafft --auto /home/heineike_wsl2/alphafold/selected_proteins/og_sequences/proteome_tm/OG1004_REF_Scer_AF-P15938-F1-model_v2.pep.fasta > /home/heineike_wsl2/alphafold/msas/sequence/mafft/fasta/OG1004_REF_Scer_AF-P15938-F1-model_v2.aln.fasta
#mafft --auto /home/heineike_wsl2/alphafold/selected_proteins/og_sequences/proteome_tm/OG1004_REF_Scer_AF-P15938-F1-model_v2.pep.fasta > /home/heineike_wsl2/alphafold/msas/sequence/mafft/fasta/OG1004_REF_Scer_AF-P15938-F1-model_v2.aln.fasta