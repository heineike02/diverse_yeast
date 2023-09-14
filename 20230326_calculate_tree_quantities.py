#import sys
import os
import subprocess
from Bio import SeqIO
import pandas as pd
import shutil
import numpy as np


base_dir = os.path.normpath('/home/heineikeb/alphafold')
aln_base_dirs = {'struct_tmalign': (base_dir + os.sep + os.path.normpath('msas/structural/tm_align'),'fasta_renamed','trees', 'tm'),
                 'seq_mafft':  (base_dir + os.sep + os.path.normpath('msas/sequence/mafft'),'fasta','trees','aln'),
                 'seq_clustalo':  (base_dir + os.sep + os.path.normpath('msas/sequence/clustalo'),'fasta','trees','aln')
                 }

#(<program>, <aln, or tree, tree_and_aln>)
phykit_programs = [('treeness', 'tree'),('saturation','tree_and_aln')]

output_file = base_dir + os.sep + os.path.normpath('msas/aln_tree_quantities.txt')

with open(output_file,'w') as f_out: 

    for aln_type, (aln_base_dir, aln_final,tree_final,first_suffix) in aln_base_dirs.items():
        print(aln_type)
        aln_dir = aln_base_dir + os.sep + aln_final
        tree_dir = aln_base_dir + os.sep + tree_final

        for (phykit_program, input_type) in phykit_programs:
            print(phykit_program)
            #Use trim_default because it contains only alignments used to make trees (in which too much trimming occurred)
            selected_alignments = [alignment_trim_fname.split('.')[0] for alignment_trim_fname in os.listdir(aln_base_dir + os.sep + 'trim_default')] 
            #selected_alignments = [fname.split('.')[0] for fname in all_tree_files if (len(fname.split('.'))==5 and (fname.split('.')[4] == 'treefile'))]
            for (jj, alignment) in enumerate(selected_alignments):
                #alignment = selected_alignments[0]
                print(str(jj) + ' ' + alignment)
                og = alignment.split('_')[0]
                
                tree_fname = tree_dir + os.sep + alignment + '.' + first_suffix + '.fasta.clipkit.treefile'

                aln_fname = aln_dir + os.sep + alignment + '.' + first_suffix + '.fasta' 

                if input_type == 'tree': 
                    phykit_cmd = ['phykit', phykit_program, tree_fname]
                elif input_type == 'aln':
                    phykit_cmd = ['phykit', phykit_program, aln_fname]
                elif input_type == 'tree_and_aln':
                    phykit_cmd = ['phykit', phykit_program, 
                                  '-a', aln_fname, '-t', tree_fname]


                #phykit_out_tmp = base_dir + os.sep + 'tmp' + os.sep + 'phykit.tmp'
                #with open(phykit_out_tmp,'w') as f_cds_trimmed:
                    #output = subprocess.run(phykit_cmd, stdout=f_cds_trimmed)
                output = subprocess.run(phykit_cmd, stdout=subprocess.PIPE)
                #Would need to process the output differently if it not just a single value
                output_val = float(output.stdout.strip())
                #mean_term_branches = float(str(output.stdout).split('\\n')[0].split(':')[1].strip())


                #print([phykit_program,output_val, og, aln_type, input_type, alignment])
                f_out.write('\t'.join([phykit_program,str(output_val), og, aln_type, input_type, alignment])+ '\n')
                
                #tree_data[alignment] = (og,mean_term_branches) #total_tree_length)
        
    #tree_data_df = pd.DataFrame.from_dict(tree_data, orient='index', columns = ['og', 'mean_term_branches']) #'total_tree_length'])
    #tree_data_df.to_csv(tree_dir + os.sep + 'tree_data.csv'
