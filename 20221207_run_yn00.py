#Running in 20221207_run_yn00.sh instead of this script

import shutil
import os
import subprocess

base_dir = os.path.normpath('/home/heineike_wsl2/alphafold')
paml_prog = '/var/lib/paml/paml-4.10.6/bin/yn00'

#Make folder to do calculation
calc_dir = base_dir + os.sep + os.path.normpath('selection_calculations/current_calc') + os.sep

og_base = 'OG4150_REF_Scer_AF-P07256-F1-model_v2'

# Copy in tree
tree_orig = base_dir + os.sep + os.path.normpath('msas/structural/tm_align/trees/' + og_base + '.tm.fasta.clipkit.treefile')
tree_renamed = tree_orig+'.renamed'
shutil.copy(tree_renamed, calc_dir + 'tree.treefile')

# Copy in phy
msa_cds_trimmed = base_dir + os.sep +  os.path.normpath('msas/structural/tm_align/cds_trim_strict/' + og_base +  '.tm.fasta.clipkit.cds')
msa_cds_trimmed_phy_codeML = msa_cds_trimmed + '.renamed.codeML.phy'
shutil.copy(msa_cds_trimmed_phy_codeML, calc_dir + 'aln.phy')

# Run codeML
subprocess.run(['/var/lib/paml/paml-4.10.6/bin/yn00', 'yn00.ctl0'])
               

# Copy out data
paml_results_out = base_dir + os.sep + os.path.normpath('selection_calculations/yn00/' + og_base + '.yn00.out')
shutil.copy(calc_dir + og_base + '.yn00.out', paml_results_out)

#clean out current directory