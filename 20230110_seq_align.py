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

base_dir = os.path.normpath('/home/heineike_wsl2/alphafold')
proteome_dir = base_dir + os.sep + os.path.normpath('selected_proteins/og_sequences/proteome_tm')

proteome_files = os.listdir(proteome_dir)
selected_proteomes = [fname.split('.')[0] for fname in proteome_files]


aln_dir = base_dir + os.sep + os.path.normpath('msas/sequence') 


#for proteome in selected_proteomes: 

proteome = selected_proteomes[0]

#Run Mafft on proteome with auto parameter.  

## Wasn't running in container. 

##Output the screen as a log so that I can extract the strategy as .aln.fasta.log
mafft_command = ['mafft', '--auto',       #'--genafpair', '--maxiterate', '1000', 
                proteome_dir + os.sep + proteome + '.pep.fasta', 
                '>',
                aln_dir + os.sep + os.path.normpath('mafft/fasta/'+ proteome + '.aln.fasta')]

print(' '.join(mafft_command))

subprocess.run(mafft_command)


#ls_cmd = ['ls', proteome_dir + os.sep]
#subprocess.run(ls_cmd)


#/home/heineike_wsl2/alphafold/selected_proteins/og_sequences/proteome_tm/OG1004_REF_Scer_AF-P15938-F1-model_v2.pep.fasta