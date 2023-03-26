#import sys
import os
import subprocess
from Bio import SeqIO
import pandas as pd
import shutil
import numpy as np


base_dir = os.path.normpath('/home/heineike_wsl2/alphafold')
tree_dir = base_dir + os.sep + os.path.normpath('msas/structural/tm_align/trees')

all_tree_files = os.listdir(tree_dir)
selected_alignments = [fname.split('.')[0] for fname in all_tree_files if (len(fname.split('.'))==5 and (fname.split('.')[4] == 'treefile'))]

phykit_programs = ['total_tree_length', 'treeness','internal_branch_stats','long_branch_score','terminal_branch_stats']

tree_data = {}



#for alignment in selected_alignments: 
alignment = selected_alignments[0]
print(alignment)
og = alignment.split('_')[0]

tree_fname = tree_dir + os.sep + alignment + '.tm.fasta.clipkit.treefile'

phykit_cmd = ['phykit', phykit_programs[1],
                  tree_fname]

phykit_out_tmp = base_dir + os.sep + 'tmp' + os.sep + 'phykit.tmp'
with open(phykit_out_tmp,'w') as f_cds_trimmed:
    output = subprocess.run(phykit_cmd, stdout=f_cds_trimmed)
#output = subprocess.run(phykit_cmd, stdout=subprocess.PIPE)

#total_tree_length = float(output.stdout.strip())
#mean_term_branches = float(str(output.stdout).split('\\n')[0].split(':')[1].strip())

#tree_data[alignment] = (og,mean_term_branches) #total_tree_length)
    
#tree_data_df = pd.DataFrame.from_dict(tree_data, orient='index', columns = ['og', 'mean_term_branches']) #'total_tree_length'])
#tree_data_df.to_csv(tree_dir + os.sep + 'tree_data.csv')



