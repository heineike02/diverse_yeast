#Using the same sequences as the structural alignment, make sequence based alignments using mafft and clustalo. 


import os
import subprocess
from Bio import SeqIO
#import pandas as pd
#import shutil
#import numpy as np

home_dir = os.path.normpath('/home/heineikeb')

base_dir = home_dir + os.sep + 'alphafold'
proteome_dir = base_dir + os.sep + os.path.normpath('selected_proteins/og_sequences/proteome_tm')

proteome_files = os.listdir(proteome_dir)
selected_proteomes = [fname.split('.')[0] for fname in proteome_files]


aln_dir = base_dir + os.sep + os.path.normpath('msas/sequence') 

#Only align those proteins that passed my filters for the sequence alignments and trees (445 in msas/structura/tm_align/cds_align

for cds_align in os.listdir(base_dir + os.sep + os.path.normpath('msas/structural/tm_align/cds_aln')):
    alignment  = cds_align.split('.')[0]
    print(alignment)

    #alignment = 'OG1254_REF_Scer_AF-P40395-F1-model_v2'  #selected_proteomes[0]

    #Run Mafft on proteome with auto parameter.  

    #Output the screen as a log so that I can extract the strategy as .aln.fasta.log
    print('Obtaining mafft alignment')
    mafft_command = ['mafft', '--auto', 
                     proteome_dir + os.sep + alignment + '.pep.fasta']

    mafft_output_fname =  aln_dir + os.sep + os.path.normpath('mafft/fasta/'+ alignment + '.aln.fasta')
    mafft_log_fname = aln_dir + os.sep + os.path.normpath('mafft/fasta/'+ alignment + '.aln.fasta.log')

    with open(mafft_output_fname, 'w') as mafft_output: 
        with open(mafft_log_fname, 'w') as mafft_log: 
            subprocess.run(mafft_command, stdout=mafft_output, stderr=mafft_log)
    

    
    #Run Clustalo on proteome with default settings
    print('Getting Clustal Alignment')
    clustalo_output_fname =  aln_dir + os.sep + os.path.normpath('clustalo/fasta/'+ alignment + '.aln.fasta')
    clustalo_log_fname = aln_dir + os.sep + os.path.normpath('clustalo/fasta/'+ alignment + '.aln.fasta.log')

    clustalo_command = [home_dir + os.sep + os.path.normpath('clustal/clustalo-1.2.4-Ubuntu-x86_64'), 
                        '-i', proteome_dir + os.sep + alignment + '.pep.fasta',
                        '-o', clustalo_output_fname,
                        '-v','-v',
                        '-l', clustalo_log_fname
                       ] 
    subprocess.run(clustalo_command)
