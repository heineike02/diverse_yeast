source('~/github/diverse_yeast/diverse_yeast_tools.R')

#For a given List of S.cer genes, make a table that includes data for all the orthologs in selected other species

#  Gal Pathway 
# Main genes upregulated in S.cer from Dalal et al 2016
#P04385, GAL1_YEAST, YBR020W
#P04397, GAL10_YEAST, YBR019C 
#P08431, GAL7_YEAST, YBR018C 
#P13181, GAL2_YEAST, YLR081W 

sc_gene_list = c('YBR020W', 'YBR019C', 'YBR018C', 'YLR081W')
list_name = 'gal_test'
source_spec = 'Scer'

other_spec_list = c('Klac', 'Lthe', 'Calb','Dhan', 'Wano', 'Spom', 'Kmar', 'Zrou', 'Kpha', 'Gcan')  # 'Ctro' - not included because of missing data

output_type = 'expNA' #either 'expNA' for raw expression data (with NA's left in), 'exp' with minimum value imputation or 'LFC'


if (output_type == 'LFC') {
    #Fold Change combinations (as calculated in proteome_analysis script)
    fc_combos = list('CN1_C2_2'= c('CN1_2', 'C2_2'),
                     'CN1_N2_2'= c('CN1_2', 'N2_2'),
                     'CN1_C2_6'= c('CN1_6', 'C2_6'),
                     'CN1_N2_6'= c('CN1_6', 'N2_6'),
                     'CN1_2_6' = c('CN1_2', 'CN1_6'),
                     'C2_2_6' = c('C2_2', 'C2_6'),
                     'N2_2_6' = c('N2_2', 'N2_6')
                     )
    

    conds = names(fc_combos)
  
  } else if (output_type %in% c('exp', 'expNA')) {
  
    conds = c('CN1_2', 'CN1_6','C2_2', 'C2_6', 'N2_2', 'N2_6')
  }

  
sc_map = data_table_from_scer_list(output_type, conds, other_spec_list)

warnings()
  



write.csv(sc_map, paste(working_dir, list_name, '.csv', sep=''))

