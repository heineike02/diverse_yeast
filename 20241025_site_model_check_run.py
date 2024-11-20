#This file checks which folders have already been run and then outputs the remaining folders
import os

base_dir = "/home/heineike/alphafold_ln"
#base_dir = os.path.normpath('G:/My Drive/Crick_LMS/projects/diverse_yeasts/alphafold')
selection_calc_dir = base_dir + os.sep + os.path.normpath('selection_calculations/site_model')

cds_list_raw = os.listdir(base_dir + os.sep + os.path.normpath('msas/structural/tm_align/cds_trim_strict'))

og_list = []

for fname in cds_list_raw: 
    fname_split = fname.split('.')
    if len(fname_split)==8:
        if (fname_split[6]=='codeML') & (fname_split[7]=='phy'):
            og_list.append(fname_split[0])

#Get list of OGs that have already been calculated

already_calc = []

for fname in os.listdir(selection_calc_dir):
    if fname[0:2]=='OG':
        already_calc.append(fname)

uncalculated = list(set(og_list)-set(already_calc))

print("There are {} uncalculated OGs.".format(len(uncalculated)))

f_out_fname = selection_calc_dir + os.sep + "uncalculated_ogs.txt"

with open(f_out_fname,'w') as f_out: 
    for og_ref in uncalculated: 
        f_out.write(og_ref + '\n')


#Also should filter out list that didn't even converge for M0
