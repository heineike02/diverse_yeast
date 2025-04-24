#This file checks which folders have already been run and then outputs the remaining folders
import os
import pandas as pd

#base_dir = "/home/heineike/alphafold_ln"
#base_dir = os.path.normpath('G:/My Drive/Crick_LMS/projects/diverse_yeasts/alphafold')
base_dir = "/home/heineikeb/alphafold"

selection_calc_dir = base_dir + os.sep + os.path.normpath('selection_calculations/site_model')

cds_list_raw = os.listdir(base_dir + os.sep + os.path.normpath('msas/structural/tm_align/cds_trim_strict'))

og_list = []

for fname in cds_list_raw: 
    fname_split = fname.split('.')
    if len(fname_split)==8:
        if (fname_split[6]=='codeML') & (fname_split[7]=='phy'):
            og_list.append(fname_split[0])

#Select only OGs that don't have dS>3 and that have no convergence issues
dnds_m0 = pd.read_csv(base_dir + os.sep + os.path.normpath('selection_calculations/m0/m0_dS_filter_20231128.csv'), index_col = 0 )

dnds_m0_alignments_present = dnds_m0.loc[list(set(dnds_m0.index) & set(og_list) ),]
#when running this initially the dS and convergence issue flag didn't work because the series were not boolean types
#dnds_m0_filt =  dnds_m0_alignments_present.loc[((dnds_m0['dS>3']==False) & (dnds_m0['convergence_issue'])==False),]  

dnds_m0_filt =  dnds_m0_alignments_present.loc[(~(dnds_m0['dS>3'].astype('bool')) & ~(dnds_m0['convergence_issue'].astype('bool'))),:]  

print("Calculating site tests for {} orthogroups without a convergence issue and with dS<=3 for M0".format(len(dnds_m0_filt)))


#Get list of OGs that have already been calculated
already_calc = []

for fname in os.listdir(selection_calc_dir):
    if fname[0:2]=='OG':
        already_calc.append(fname)

uncalculated = list(set(dnds_m0_filt.index)-set(already_calc))
#uncalculated = list(set(og_list)-set(already_calc))

print("There are {} uncalculated OGs.".format(len(uncalculated)))

f_out_fname = selection_calc_dir + os.sep + "uncalculated_ogs.txt"

with open(f_out_fname,'w') as f_out: 
    for og_ref in uncalculated: 
        f_out.write(og_ref + '\n')


#Also should filter out list that didn't even converge for M0
