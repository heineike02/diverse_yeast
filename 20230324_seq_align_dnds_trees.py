#Make trees based on default trimmed alignments from sequence-based MSAs generated with MAFFT and CLUSTALO

#import sys
import os
import subprocess
from Bio import SeqIO
import pandas as pd
import shutil

base_dir = os.path.normpath('/home/heineikeb/alphafold')
#os.path.normpath('/home/heineike_wsl2/alphafold') #Ben's computer

for aln_type in ['mafft', 'clustalo']:

    aln_dir = base_dir + os.sep + os.path.normpath('msas/sequence/' + aln_type) 

    align_files = os.listdir(aln_dir + os.sep + 'fasta')  
    selected_alignments_all = [fname.split('.')[0] for fname in align_files]

    completed_alignments = [fname.split('.')[0] for fname in os.listdir(aln_dir + os.sep + 'trim_default')]

    selected_alignments = list(set(selected_alignments_all)-set(completed_alignments))

    #os.path.normpath('/home/heineike_wsl2/Crick_LMS/projects/diverse_yeasts/alphafold')
    # output_dir = base_dir + os.sep + os.path.normpath('selection_calculations/20220526_sel_calc')

    #selected_ogs = ['OG2645'] #'OG4150', 'OG2603', 'OG3677', 'OG2845']

    #selected_og_refs = ['OG4150_REF_Scer_AF-P07256-F1-model_v2'] #, 'OG2603_REF_Scer_AF-P50076-F1-model_v2', 'OG2845_REF_Scer_AF-P43577-F1-model_v2', 'OG3677_REF_Scer_AF-P47125-F1-model_v2', 'OG1299_REF_Scer_AF-P00549-F1-model_v2']

    min_seq = 4

    trees_log_fname = aln_dir + os.sep + os.path.normpath('trees/trees_log.txt')

    with open(trees_log_fname, 'w') as trees_log:

        for jj, alignment in enumerate(selected_alignments): 
            print(alignment + ' ' + str(jj) + ' of ' + str(len(selected_alignments)))
            #og_ref = 'OG4150_REF_Scer_AF-P07256-F1-model_v2'
            #og,ref = og_ref.split('_REF_')
            
            print('Trimming Alignment with default parameters')
            og = alignment.split('_')
            og_pep_msa_fname = aln_dir + os.sep + 'fasta' + os.sep + alignment + '.aln.fasta'

            #Verify that there are only sequences greater than or equal to min_seq. 
            og_pep_msa = SeqIO.parse(og_pep_msa_fname,'fasta')
            nseqs = len(list(og_pep_msa))    
            
            assert nseqs>=min_seq, alignment + ' has less than ' + str(min_seq) + 'sequences' 

            #trim alignment with default settings and outputs log file
            clipkit_cmd = ['clipkit', og_pep_msa_fname, '-l']
            subprocess.run(clipkit_cmd)

            #Move log and output

            suffixes = ['.clipkit','.clipkit.log' ]

            for suffix in suffixes: 
                #Move clipkit files to new folder
                shutil.move(og_pep_msa_fname + suffix, 
                                aln_dir + os.sep + os.path.normpath('trim_default/' + alignment + '.aln.fasta' + suffix)
                               )

                og_pep_msa_fname_trimmed = aln_dir + os.sep + os.path.normpath('trim_default/' + alignment + '.aln.fasta.clipkit')

            # Run iQtree on trimmed peptide MSA  
            # Should I run with a pombe outgroup? 
            
            print('Building Protein Tree')

            iqtree_command = ["iqtree", 
                              "-s" , og_pep_msa_fname_trimmed,
                              #"-m", 'LG+I+G4',  #'MF', #only runs model finder 
                              "-nt", "7", #"AUTO"  automatically determines number of threads but 7 was performing well
                              "-bb", "1000",
                              "-alrt", "1000",
                              #"-o", 'Spom_AF-Q10208-F1-model_v2'  #Outgroup for rooting should be pombe  for now using default. 
                             ]
            #print(" ".join(iqtree_command))

            subprocess.run(iqtree_command)

            #move treefiles to new directory
            
            tree_output_files = os.listdir(aln_dir + os.sep + os.path.normpath('trim_default'))

            for suffix in [ 'ckp.gz','iqtree', 'bionj','mldist', 'log', 'treefile', 'contree','model.gz','splits.nex','uniqueseq.phy']:
                fname_from = alignment + '.aln.fasta.clipkit.' + suffix

                if fname_from in tree_output_files: 
                    fname_from_full = aln_dir + os.sep + os.path.normpath('trim_default/' + fname_from)            
                    fname_to = aln_dir + os.sep + os.path.normpath('trees/' + alignment + '.aln.fasta.clipkit.' + suffix)
                    shutil.move(fname_from_full, fname_to)
                
            #Format phylogenetic Tree for codeml by shortening the name

            # shorten name.  Uses seq_name_map made for tm-align alignments in 20221206_struct_align_dnds_msas.py
            print('Formatting tree for Codeml')
            seq_name_map_fname = aln_dir + os.sep +  os.path.normpath('seq_name_map/' + alignment + '.tm.tsv')

            #Use phykit rename_tree_tips to shorten the name
            tree_orig = aln_dir + os.sep + os.path.normpath('trees/' + alignment + '.aln.fasta.clipkit.treefile')
            tree_renamed = tree_orig+'.renamed'
            phykit_rename_cmd = ['phykit', 'rename_tree_tips',
                             tree_orig,
                             '-i', seq_name_map_fname, 
                             '-o', tree_renamed
                            ]

            #print(" ".join(phykit_rename_cmd))

            subprocess.run(phykit_rename_cmd)

            #else: 
            #   trees_log.write(alignment + ' has less than ' + str(min_seq) + ' sequences. No tree created.\n')

