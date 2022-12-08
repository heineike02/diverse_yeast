#import json
import sys
import os
import subprocess
from Bio import SeqIO
import pandas as pd

## To install clipkit, phykit and biokit, first created a virtual environment per the following commands below: 

## Create virtual environment
# python -m venv .venv_biokit
## activated virtual environment
# source .venv/bin/activate
## Installed software
# pip install clipkit
# pip install phykit
# pip install biokit

#begin by activating source: 
#source /home/heineike/.venv_biokit/bin/activate

base_dir = os.path.normpath('/home/heineike/Crick_LMS/projects/diverse_yeasts/alphafold')
# output_dir = base_dir + os.sep + os.path.normpath('selection_calculations/20220526_sel_calc')

#clipkit <input_file> -l 
#G:\My Drive\Crick_LMS\projects\diverse_yeasts\alphafold\msas\structural\fasta_filt
#trim alignment and output log file

og_pep_fname =  'OG1254_REF_Scer_AF-P40395-F1-model_v2.struct_filt.fasta'   #'OG1299_REF_Scer_AF-P00549-F1-model_v2.struct_filt.fasta'
og_ref_base = og_pep_fname.split('.')[0]
#og,ref = og_ref_base.split('_REF_')

msa_pep_fname = base_dir + os.sep + os.path.normpath('msas/structural/fasta_filt/' + og_pep_fname)
clipkit_cmd = ['clipkit', msa_pep_fname, '-l']

#Toni recommended this command from trimal: ./trimal -in OG1299_cds_aln.pro -backtrans OG1299_cds_aln.fst
#-strictplus -phylip_paml   -out OG1299_cds_aln_trimmed.cds

subprocess.run(clipkit_cmd)

#Move the clipkit log and trimmed alignment to a new place

#cds_fname = base_dir + os.sep +  os.path.normpath('og_sequences/cds/' + og_ref_base + '.cds.fasta')
msa_cds_fname = base_dir + os.sep +  os.path.normpath('msas/structural/fasta_filt_cds/' + og_ref_base + '.struct_filt_cds.fasta')
msa_cds_trimmed = base_dir + os.sep +  os.path.normpath('msas/structural/fasta_filt_cds_trimmed/' + og_ref_base + '.struct_filt_cds.trimmed.fasta')

#Note:  For phykit need to provide corresponding alignments and the clipkit log for the trimmed alignment.  
phykit_cmd = ['phykit', 'thread_dna',
             '-p', msa_pep_fname + '.clipkit',
             '-n', msa_cds_fname, 
             '-c', msa_pep_fname + '.clipkit.log'
            ]

with open(msa_cds_trimmed,'w') as f_cds_trimmed:
    output = subprocess.run(phykit_cmd, stdout=f_cds_trimmed)
#phykit thread_dna -p <file> -n <file> -l <clipkit_log>




# #Use IQTree to make tree for the trimmed file
# iqtree_command = ["/home/heineike/iqtree-1.6.12-Linux/bin/iqtree", 
#                   "-s" , msa_pep_fname,
#                   #"-m", 'MF', #only runs model finder 
#                   "-nt", "AUTO",  #automatically determines number of threads 
#                   "-o", 'Spom_AF-Q10208-F1-model_v2']
# print(" ".join(iqtree_command))
#Model LG+I+G4 was the best by BIC, LG+R5 did well by AIC and corrected AIC.  Use LG+I+G4

# seq_name_map_fname = base_dir + os.sep +  os.path.normpath('og_sequences/seq_name_map/' + og_ref_base + '.tsv')

# tree_renamed_out = base_dir + os.sep +  os.path.normpath('msas/structural/tree_renamed/' + og_ref_base + '.clipkit.renamed.treefile')
# phykit_rename_cmd = ['phykit', 'rename_tree_tips',
#              msa_pep_fname + '.clipkit.treefile',
#              '-i', seq_name_map_fname, 
#              '-o', tree_renamed_out
#             ]

# subprocess.run(phykit_rename_cmd)



# def reformat_tree(tree):
#     alle=open(tree).read()
#     trees=[]
#     for a in range(len(alle)):
#         if a==0: continue
#         if alle[a-1]=='(' and alle[a]!='(':
#             left,right=alle[a:].split(':',1)
#             trees.append(alle[:a]+left+" #1:"+right)

#         if alle[a-1]==')':
#             left,right=alle[a-1],alle[a-1:]
#             trees.append(left+" #1:"+right)

#         if a<len(alle):
#             if alle[a]==',' and alle[a+1]!='(':
#                 left,right=alle[a+1:].split(':',1)
#                 #print(name)
#                 trees.append(alle[:a+1]+left+" #1:"+right)
#     return(trees)





# newtrees=reformat_tree(tree_renamed_out)
# tree_reformatted_out = base_dir + os.sep +  os.path.normpath('msas/structural/tree_renamed/' + og_ref_base + '.clipkit.renamed.reform.treefile')

# with open(tree_reformatted_out,'w') as f_out_tree_reform:
#     for trees in range(len(newtrees)):
#         f_out_tree_reform.write(newtrees[trees])
    






# # #Rename cds alignment 

# msa_cds_trimmed_renamed = base_dir + os.sep +  os.path.normpath('msas/structural/fasta_filt_cds_trimmed/' + og_ref_base + '.struct_filt_cds.trimmed.renamed.fasta')

# #shorten the name of the fasta and build name_replace dictionary
# seq_name_map_df = pd.read_csv(seq_name_map_fname, sep='\t')
# seq_name_map = dict(zip(seq_name_map_df['seq_name'],seq_name_map_df['seq_no']))

# msa_cds_trimmed_seqs = SeqIO.parse(msa_cds_trimmed, 'fasta')
# with open(msa_cds_trimmed_renamed,'w') as f_out: 
#     for record in msa_cds_trimmed_seqs:
#         new_id = seq_name_map[record.id] 
#         f_out.write('>' + str(new_id) + '\n')
#         f_out.write(str(record.seq) + '\n')


# #converts fasta to PHYLIP format 
# msa_cds_trimmed_phy = base_dir + os.sep +  os.path.normpath('msas/structural/fasta_filt_cds_trimmed/' + og_ref_base + '.struct_filt_cds.trimmed.renamed.phy')
# biokit_cmd = ['biokit','file_format_converter',
#               '-i',msa_cds_trimmed_renamed, '-iff', 'fasta',
#               '-o', msa_cds_trimmed_phy,
#               '-off', 'phylip']
              
# subprocess.run(biokit_cmd)

# #adds I in phylip to indicate interleaved status
# msa_cds_trimmed_phy_codeML = base_dir + os.sep +  os.path.normpath('msas/structural/fasta_filt_cds_trimmed/' + og_ref_base + '.struct_filt_cds.trimmed.renamed.codeML.phy')

# with open(msa_cds_trimmed_phy,'r') as f_in:
#     with open(msa_cds_trimmed_phy_codeML,'w') as f_out: 
#         line = next(f_in)
#         line_out = line.strip('\n') + ' I\n'
#         f_out.write(line_out)
#         for line in f_in:
#             f_out.write(line)




# build codeML file



# run codeML








#Input: list of OGs, (?parameters for codeML), specific protein ids to remove (as .json), suffix for calculation (defoult '')
# #og_input = {'OG1299': {'name': 'CDC19_PYK2',
#                           'genes_to_rm': ['0_2015']
#                          },
#                'OG1355': {'name': 'ERG11',
#                           'genes_to_rm': ['0_2015']
#                          },
#                'OG1390': {'name': 'STR2_HSU1',
#                           'genes_to_rm': ['0_2015']
#                          }
#               }

#Output: CodeML calculations for each OG

##eventually use sys.argv to have command line input

# #Load input file
# base_dir = os.path.normpath('/home/heineike/Crick_LMS/projects/diverse_yeasts/alphafold')
# output_dir = base_dir + os.sep + os.path.normpath('selection_calculations/20220526_sel_calc')


# if not(os.path.isdir(output_dir)):
#     os.mkdir(output_dir)
# else: 
#     print('Directory ' + output_dir + ' already exists')

# with open(base_dir + os.sep + os.path.normpath('selection_calculations/20220526_og_input.json'), 'r') as og_input_f:
#     og_input  = json.load(og_input_f)

# #For each OG


# #for og, param_dict in og_input.items(): 
# og = 'OG1299'
# genes_to_rm = og_input[og]['genes_to_rm']

# #makes directory
# og_dir = output_dir + os.sep + og
# if not(os.path.isdir(og_dir)):
#     os.mkdir(og_dir)
# else: 
#     print('Directory ' + og_dir + ' already exists')




# #Load Tree



# with open(base_dir + os.sep + os.path.normpath('alphafold/selection_calculations/20220526_og_input.json'), 'w') as og_input_f: 
#     og_input = json.dump(og_input, fp = og_input_f)

#fix names
#Saves in new directory

#writes codeml.ctl file

#runs codeml