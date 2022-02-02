source('~/github/diverse_yeast/diverse_yeast_tools.R')


make_pairplots = FALSE   #Flag to produce pairplots for each species
save_species_data = TRUE #Flag to save 

specs = c('Scer-BY4741KI', 'Scer', 'Klac', 'Lthe', 'Calb', 'Ctro', 'Dhan', 'Wano', 'Spom', 'Kmar', 'Zrou', 'Kpha', 'Gcan')

#This might not be needed
spec4_to_spec2 = c('Scer-BY4741KI'='BY', 'Scer'='SC', 'Klac'='KL', 'Lthe'='LT', 'Calb'='CA', 'Ctro'='CT', 'Dhan'='DH', 'Wano'='WA', 'Spom'='SP', 'Kmar'='KM', 'Zrou'='ZR', 'Kha'='PP', 'Gcan'='GC')
#The species name for the uniprot file we used for picchia pastoris is Komatagella phaffii, so I think we should use that for the 3 letter code.

#For the s. cerevisiae data it was saved as an R. File.  Opened it and then Saved it as a .tsv using
#write.table(proteinWide_woQC, paste(working_dir, spec, '/','BY_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_ProteinIds_woQC.tsv', sep=''))

#Need to explicitly write the name of the source data file
spec_data_fnames = c('Scer-BY4741KI' = 'BY_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_woQC.tsv', 
                     'Scer' = 'SC_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_woQC.tsv',
                     'Klac' = 'KL_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_woQC.tsv',
                     'Lthe' = 'LT_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_woQC.tsv',
                     'Calb' = 'CA_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_woQC.tsv',
                     'Ctro' = 'CT_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_woQC.tsv',
                     'Dhan' = 'DH_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_woQC.tsv',
                     'Wano' = 'WA_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_woQC.tsv',
                     'Spom' = 'SP_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_woQC.tsv',
                     'Kmar' = 'KM_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_woQC.tsv',
                     'Zrou' = 'ZR_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_woQC.tsv',
                     'Kpha' = 'PP_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_woQC.tsv',
                     'Gcan' = 'GC_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_woQC.tsv'
)



narep_pct_min = 0.9
infrep_pct_max = 1.05

cond_order = c('CN1_2', 'CN1_6','C2_2', 'C2_6', 'N2_2', 'N2_6')

#Fold Change is calculated for various pairs of gene expression value
fc_combos = list('CN1_C2_2'= c('CN1_2', 'C2_2'),
                 'CN1_N2_2'= c('CN1_2', 'N2_2'),
                 'CN1_C2_6'= c('CN1_6', 'C2_6'),
                 'CN1_N2_6'= c('CN1_6', 'N2_6'),
                 'CN1_2_6' = c('CN1_2', 'CN1_6'),
                 'C2_2_6' = c('C2_2', 'C2_6'),
                 'N2_2_6' = c('N2_2', 'N2_6')
)



#For each species saves fold change data and visualizes matrix of scatter plots ussing ggpairs
exp_list_na = list()
exp_list = list()
fc_list = list()

for (spec in specs) {
  cat(paste(spec), '\n')
  #Load protein table for each species
  protein_data_fname = paste(working_dir, spec, '/', spec_data_fnames[[spec]], sep='')
  protein_data = read.table(file=protein_data_fname, header=TRUE)
  rownames(protein_data) = protein_data$Protein.Group
  protein_data$Protein.Group = NULL
  
  #Rename columns to remove species and sample info
  new_cols = sapply(colnames(protein_data), col_xform)
  colnames(protein_data) = new_cols
  protein_data = protein_data[, intersect(cond_order, colnames(protein_data))]
  
  #save data with NA (without minimum value imputation)
  protein_data_log_NA = log10(protein_data)
  
  exp_list_na[[spec]] = protein_data_log_NA
  
  
  
  #NA Replacement Rule: 
  #Remove rows that have NAs for all conditions
  narow_count = rowSums(as.data.frame(lapply(as.data.frame(is.na(protein_data)), as.numeric)))
  if (length(which(narow_count == length(colnames(protein_data)))) > 0){ 
    print('Row of Data with all NAs - need to add code to remove')
  }
  #Replace NAs in other rows with minimum value of 90% of the minimum value in the rest of the data set. 
  narep = narep_pct_min*min(protein_data, na.rm=TRUE)
  narep_by_col = setNames(as.list(rep(narep, length(colnames(protein_data)))), colnames(protein_data))
  protein_data = tidyr::replace_na(protein_data, narep_by_col ) #needed to add tidyr::
  
  #Replace Inf with 105% of the maximum value in the rest of the data set
  infrep = infrep_pct_max * max(as.matrix(protein_data[is.finite(as.matrix(protein_data))]))
  protein_data[is.infinite(as.matrix(protein_data))]=infrep
  
  
  
  protein_data_log = log10(protein_data)
  

  
  exp_list[[spec]] = protein_data_log
  
  #make ggpairs plot
  if (make_pairplots) {
    plotrange = c(2.5,7)
    p_exp = ggpairs(protein_data_log,
                  lower = list(continuous = wrap(lowerfun, pt_alpha=0.2, pt_size=1, plotrange=plotrange)),
                  diag = list(continuous = wrap(diagfun, plotrange = plotrange ))
  )
    p_exp = p_exp + labs(title= paste(spec, 'Log Expression'))
    show(p_exp)
  }
  
  #Calculate LFC
  
  protein_data_fc = data.frame(row.names=rownames(protein_data_log))
  for (fc_combo in names(fc_combos)) {
    if (!(length(intersect(fc_combos[[fc_combo]], colnames(protein_data_log)))==2)) {
      print(paste('One or both conditions from ', fc_combo, ' not present in data.  Species: ', spec))
    } else {
      protein_data_fc[fc_combo]=protein_data_log[fc_combos[[fc_combo]][2]]-protein_data_log[fc_combos[[fc_combo]][1]]
    }
  }
  
  fc_list[[spec]] = protein_data_fc
  
  #make ggpairs plot of different LFC comparisons
  if (make_pairplots) {
    plotrange = c(-2,2)
    p_fc = ggpairs(protein_data_fc,
                 lower = list(continuous = wrap(lowerfun, pt_alpha=0.2, pt_size=1, plotrange=plotrange)),
                 diag = list(continuous = wrap(diagfun, plotrange = plotrange ))
  )
  
    p_fc = p_fc + labs(title=paste(spec,'Fold Change') )
  
    show(p_fc)
  }
  
}

#Save expression and LFC data as .csv files
if (save_species_data) {
  for (spec in specs) {
    write.csv(exp_list[[spec]], paste(working_dir, spec, '/exp_data_',spec,'.csv', sep=''))
    write.csv(exp_list_na[[spec]], paste(working_dir, spec, '/expNA_data_',spec,'.csv', sep=''))
    write.csv(fc_list[[spec]], paste(working_dir, spec, '/LFC_data_',spec,'.csv', sep=''))
  }
}

