import json
import sys
import os
import subprocess

#begin by activating source: 
#source /home/heineike/.venv_biokit/bin/activate

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

#Load input file
base_dir = os.path.normpath('/home/heineike/Crick_LMS/projects/diverse_yeasts/alphafold')
output_dir = base_dir + os.sep + os.path.normpath('selection_calculations/20220526_sel_calc')


if not(os.path.isdir(output_dir)):
    os.mkdir(output_dir)
else: 
    print('Directory ' + output_dir + ' already exists')

with open(base_dir + os.sep + os.path.normpath('selection_calculations/20220526_og_input.json'), 'r') as og_input_f:
    og_input  = json.load(og_input_f)

#For each OG


#for og, param_dict in og_input.items(): 
og = 'OG1299'
genes_to_rm = og_input[og]['genes_to_rm']

#makes directory
og_dir = output_dir + os.sep + og
if not(os.path.isdir(og_dir)):
    os.mkdir(og_dir)
else: 
    print('Directory ' + og_dir + ' already exists')

#Loads MSA
#shortens names
#removes filter genes
#Saves in new directory

#shorten the name of the fasta and build name_replace dictionary
fasta_in = base_dir + os.sep + os.path.normpath('msas/FILES_ogs_cds_threaded/' + og + '.mfna.threaded')
fasta_short = og_dir + os.sep + og + '_aln.fasta'
name_replace = {}
with open(fasta_in,'r') as f_in:
    with open(fasta_short,'w') as f_out: 
        for line in f_in: 
            if line[0]=='>':
                new_id = line.split('|')[1]
                line_out = '>' + new_id
                name_replace[line] = new_id
            else: 
                line_out = line
            f_out.write(line_out)

fasta_filt = og_dir+os.sep+og+ '_aln_filt.fasta'
with open(fasta_short,'r') as f_in:
    with open(fasta_filt,'w') as f_out: 
        for line in f_in:
            if (line[0] == '>'):
                gene_id = line.split('>')[1].strip()
                #print(genes_to_rm)
                #print(gene_id in genes_to_rm)
                if gene_id in genes_to_rm:
                    next(f_in)  #skips the sequence line as well
                else: 
                    f_out.write(line)
            else: 
                f_out.write(line)

#converts fasta to PHYLIP format 
phy = og_dir + os.sep + og + '_aln.phy'
biokit_cmd = ['biokit','file_format_converter',
              '-i',fasta_filt, '-iff', 'fasta',
              '-o', phy,
              '-off', 'phylip']
              
subprocess.run(biokit_cmd)

#adds I in phylip to indicate interleaved status
phy = og_dir + os.sep + og + '_aln.phy'
phy_codeML = og_dir + os.sep + og + '_codeML_aln.phy'

with open(phy,'r') as f_in:
    with open(phy_codeML,'w') as f_out: 
        line = next(f_in)
        line_out = line + ' I'
        f_out.write(line_out)
        for line in f_in:
            f_out.write(line_out)


#Load Tree



with open(base_dir + os.sep + os.path.normpath('alphafold/selection_calculations/20220526_og_input.json'), 'w') as og_input_f: 
    og_input = json.dump(og_input, fp = og_input_f)

#fix names
#Saves in new directory

#writes codeml.ctl file

#runs codeml