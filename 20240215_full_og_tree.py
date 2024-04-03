#Make trimmed alignments in order to display a phylogenetic tree for an orthogroup before it was split into subgroups.  
#
#First make a standard trimming with clipkit
#
#Then run IQtree 
#
#Input: og
#
#Output: 
#  Strictly trimmed protein alignment and codon alignment (used as input for codemL)
#  Phylip file used as input for Codon alignment. 

#import sys
import os
import subprocess
from Bio import SeqIO
import pandas as pd
import shutil

#base_dir = os.path.normpath('/home/heineikeb/alphafold') #MS03 server
base_dir = os.path.normpath('/home/heineike/alphafold') #Ben's computer

aln_dir =  base_dir + os.sep + os.path.normpath('selected_proteins/pdbs/OG1022') 

#use /tm_align/cds_aln to pick the aligments as that is post filtering for alignments that become very short upon trimming
aln_file = aln_dir + os.sep + 'us_align_clean.fasta' 

#tree_log_fname = aln_dir + os.sep + os.path.normpath('tree/tree_log.txt')
#with open(tree_log_fname, 'w') as trees_log:

print('Trimming and tree creation for ' + aln_file)
#og_ref = 'OG4150_REF_Scer_AF-P07256-F1-model_v2'
#og,ref = og_ref.split('_REF_')

print('Trimming Alignment with default parameters')

#trim alignment with default settings and outputs log file
clipkit_cmd = ['clipkit', aln_file, '-l']
subprocess.run(clipkit_cmd)

#Move log and output

suffixes = ['.clipkit','.clipkit.log' ]

for suffix in suffixes: 
    #Move clipkit files to new folder
    shutil.move(aln_file + suffix, 
                aln_dir + os.sep + os.path.normpath('tree/us_align_trimmed.fasta' + suffix)
                )

aln_fname_trimmed = aln_dir + os.sep + os.path.normpath('tree/us_align_trimmed.fasta.clipkit')

# Run iQtree on trimmed peptide MSA  
# Should I run with a pombe outgroup? 

print('Building Protein Tree')

iqtree_command = ["iqtree", 
                  "-s" , aln_fname_trimmed,
                  #"-m", 'LG+I+G4',  #'MF', #only runs model finder 
                  "-nt", "7", #"AUTO"  automatically determines number of threads but 7 was performing well
                  "-bb", "1000",
                  "-alrt", "1000",
                  #"-o", 'Spom_AF-Q10208-F1-model_v2'  #Outgroup for rooting should be pombe  for now using default. 
                 ]
#print(" ".join(iqtree_command))

subprocess.run(iqtree_command)

#move treefiles to new directory

#tree_output_files = os.listdir(aln_dir + os.sep + os.path.normpath('trim_default'))

#for suffix in [ 'ckp.gz','iqtree', 'bionj','mldist', 'log', 'treefile', 'contree','model.gz','splits.nex','uniqueseq.phy']:
#    fname_from = alignment + '.tm.fasta.clipkit.' + suffix
#
#    if fname_from in tree_output_files: 
#        fname_from_full = aln_dir + os.sep + os.path.normpath('trim_default/' + fname_from)            
#        fname_to = aln_dir + os.sep + os.path.normpath('trees/' + alignment + '.tm.fasta.clipkit.' + suffix)
#        shutil.move(fname_from_full, fname_to)
    
