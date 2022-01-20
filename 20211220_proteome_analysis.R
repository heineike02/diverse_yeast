library(ggplot2)
library(GGally)
library(plotly)

working_dir = "/camp/home/heineib/working/Ben/diverse_strains/processed_data/"

specs = c('Scer', 'Klac', 'Zrou')

#This might not be needed
spec4_to_spec2 = c('Scer'='SC','Klac'='KL', 'Zrou'='ZR')

#For the s. cerevisiae data it was saved as an R. File.  Opened it and then Saved it as a .tsv using
#write.table(proteinWide_woQC, paste(working_dir, spec, '/','BY_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_ProteinIds_woQC.tsv', sep=''))

#Need to explicitly write the name of the source data file
spec_data_fnames = c('Scer' = 'BY_ProteinWide_BatchCorrected_SF0MinPrecNum3_Stringent_ProteinIds_woQC.tsv', 
                     'Klac' = 'KL_ProteinWide_BatchCorrected_0_CV0_Stringent_woQC.tsv',
                     'Zrou' = 'ZR_ProteinWide_BatchCorrected_0_CV0_Stringent_woQC.tsv'
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


#Function for lower scatter plots to keep plot range the same for X and Y for all plots. 
lowerfun <- function(data,mapping, pt_alpha, pt_size, plotrange){
  ggplot(data = data, mapping = mapping)+
  geom_point(alpha=pt_alpha, size=pt_size) +
  scale_x_continuous(limits = plotrange)+
  scale_y_continuous(limits = plotrange)
}

#Function to pass plotrange parameters to distribution plots on the center diagonal.  
diagfun <- function(data, mapping, plotrange){
  ggplot(data = data, mapping = mapping) +
  geom_density(data=data, mapping=mapping) +
  xlim(plotrange)
}

#Extracts condition column names from long names (that include species name)
col_xform = function(col) {
  colsplit = strsplit(col, '_')
  new_col = paste(colsplit[[1]][4:5], collapse='_')
}



#For each species saves fold change data and visualizes matrix of scatter plots ussing ggpairs
exp_list = list()
fc_list = list()

for (spec in specs) {
  #Load protein table for each species
  protein_data_fname = paste(working_dir, spec, '/', spec_data_fnames[[spec]], sep='')
  protein_data = read.table(file=protein_data_fname, header=TRUE)
  rownames(protein_data) = protein_data$Protein.Group
  protein_data$Protein.Group = NULL
  
  #NA Replacement Rule: 
  #Remove rows that have NAs for all conditions
  narow_count = rowSums(as.data.frame(lapply(as.data.frame(is.na(protein_data)), as.numeric)))
  if (length(which(narow_count == length(colnames(protein_data)))) > 0){ 
    print('Row of Data with all NAs - need to add code to remove')
  }
  #Replace NAs in other rows with minimum value of 90% of the minimum value in the rest of the data set. 
  narep = narep_pct_min*min(protein_data, na.rm=TRUE)
  narep_by_col = setNames(as.list(rep(narep, length(colnames(protein_data)))), colnames(protein_data))
  protein_data = replace_na(protein_data, narep_by_col )
  
  #Replace Inf with 105% of the maximum value in the rest of the data set
  infrep = infrep_pct_max * max(as.matrix(protein_data[is.finite(as.matrix(protein_data))]))
  protein_data[is.infinite(as.matrix(protein_data))]=infrep
  
  
  
  protein_data_log = log10(protein_data)
  
  #Rename columns to remove species and sample info
  new_cols = sapply(colnames(protein_data_log), col_xform)
  colnames(protein_data_log) = new_cols
  protein_data_log = protein_data_log[, cond_order]
  
  exp_list[[spec]] = protein_data_log
  
  #make ggpairs plot
  plotrange = c(2.5,7)
  p_exp = ggpairs(protein_data_log,
                  lower = list(continuous = wrap(lowerfun, pt_alpha=0.2, pt_size=1, plotrange=plotrange)),
                  diag = list(continuous = wrap(diagfun, plotrange = plotrange ))
                  )
  
  p_exp = p_exp + labs(title= paste(spec, 'Log Expression'))
  show(p_exp)
  
  #Calculate LFC
  
  protein_data_fc = data.frame(row.names=rownames(protein_data_log))
  for (fc_combo in names(fc_combos)) {
    protein_data_fc[fc_combo]=protein_data_log[fc_combos[[fc_combo]][2]]-protein_data_log[fc_combos[[fc_combo]][1]]
  }
  
  fc_list[[spec]] = protein_data_fc
  
  #make ggpairs plot of different LFC comparisons
  plotrange = c(-2,2)
  p_fc = ggpairs(protein_data_fc,
                lower = list(continuous = wrap(lowerfun, pt_alpha=0.2, pt_size=1, plotrange=plotrange)),
                diag = list(continuous = wrap(diagfun, plotrange = plotrange ))
                )
  
  p_fc = p_fc + labs(title=paste(spec,'Fold Change') )
  
  show(p_fc)

  }



#Calculate LFC


#make ggpairs plot of different LFC comparisons


#Compare different LFC comparisons across species using ortholog mapping

#Start with 'CN1_C2_2'

output_cond = 'CN1_C2_2'
specA = 'Scer'
specB = 'Zrou'
output_type = 'LFC'
output_list = list('LFC'= fc_list, 'exp'=exp_list)

outputA_all = output_list[[output_type]][[specA]] 
outputB_all = output_list[[output_type]][[specB]] 

#Load Ortholog Mapping
orth_map = read.csv(file=paste(working_dir, 'ortholog_maps/', specB, '_', specA, '.csv', sep=''), header=TRUE, row.names=1 )  

orth_map_source_genename = strsplit(orth_map$source_genename, '[|]')


orth_map_name_parse = function(name_split) {
  name_out = strsplit(name_split[3], '_')[[1]][1]
  return(name_out)
}

orth_map_source_new_names = sapply(orth_map_source_genename, orth_map_name_parse)

orth_map$source_genename_short = orth_map_source_new_names

#Check all output genes are included in ortholog map
if (nrow(outputB_all) != length(intersect(rownames(outputB_all), orth_map$source_genename_short))) {
  warning('Some outputs not present in ortholog map')
}


#Load Protein Annotation for S. cerevisiae and use it to map orf names onto S. cer data

#saved annotations_proteins as a .tsv using data from 210318_annotation_tables.Rdata
#write.table(annotations_proteins, paste(working_dir,'annotations_proteins.tsv', sep = ''))

scer_annotation = read.delim(file=paste(working_dir, 'uniprot-proteome_UP000002311.tab', sep=''),sep='\t', header=TRUE)

#Switched to uniprot annotation file because it was missing some assignments from uniprot gene name to uniprot ID
#scer_annotation = read.table(file=paste(working_dir, 'annotations_proteins.tsv', sep=''), header=TRUE)

#renames columns to match original scer_annotation columns
scer_annotation = rename(scer_annotation, c('uniprot_id'='Entry', 'systematic_name'='Gene.names...ordered.locus..'))

print("List of Uniprot IDs present in the data that are not in the annotation file:")
print(setdiff(rownames(fc_list[['Scer']]), scer_annotation$uniprot_id))

rownames(scer_annotation) = scer_annotation$systematic_name


#Code for old annotation_proteins.tsv file.  
#
#rownames(scer_annotation) = scer_annotation$systematic_name


#Assess items in the annotation file that have duplicate uniprot ids and are NA for uniprot ID, but have a genename assigned.
uniprot_split = strsplit(scer_annotation$uniprot_id, ',')
uniprot_lengths = lapply(uniprot_split,length)
scer_annotation_uniprot_doubles = scer_annotation[which(uniprot_lengths==2),]
scer_annotation_uniprot_na = scer_annotation[which(is.na(scer_annotation$uniprot_id)),]
scer_annotation_uniprot_na_genename_present = scer_annotation_uniprot_na[which(!is.na(scer_annotation_uniprot_na$gene_name)),]
write.table(scer_annotation_uniprot_doubles, paste(working_dir, 'scer_duplicate_uniprot_ids.tsv', sep=''))
write.table(scer_annotation_uniprot_na_genename_present, paste(working_dir, 'scer_uniprot_na_genename_present.tsv', sep=''))


#The New Uniprot id mapping doesn't have the double uniprot names listed and doesn't have NAs for Uniprot when a genename is present. 

#There's no way we could recover data for the genes that have NA for Uniprot but data linked to genename - I assume DIANN throws those out
#Not sure how diann deals with duplicate uniprot ids



output_comb = data.frame()
#colnames(orth_map_comb) = c('genename_B', 'LFC_B', 'genename_A', 'sc_orf', 'LFC_A', 'orth_type')

for (orth_type in c('no_eggnog_orthologs','no_target_orthologs' )) {
  orth_map_type_subset = orth_map[which(orth_map$orth_type==orth_type),]
  output_comb_type = data.frame(genename_B = orth_map_type_subset$source_genename_short, 
                                LFC_B = outputB_all[orth_map_type_subset$source_genename_short, output_cond],
                                genename_A = 'NONE', 
                                sc_orf = 'NONE', 
                                LFC_A = 0, 
                                orth_type = orth_type
                                )
  
  output_comb = rbind(output_comb, output_comb_type)
}

for (orth_type in c('one2one', 'no_target_orthologs', 'one2many', 'many2one', 'many2many')) {
  orth_map_type_subset = orth_map[which(orth_map$orth_type==orth_type),]
  output_comb_type = data.frame(genename_B = orth_map_type_subset$source_genename_short, 
                                LFC_B = outputB_all[orth_map_type_subset$source_genename_short, output_cond],
                                orth_type = orth_type, 
                                sc_orf = orth_map_type_subset$target_genename, 
                                LFC_A = NA
  )
  
  #output_comb_type$genename_A = scer_annotation[orth_map_type_subset$target_genename, c('gene_name')]  #If data is mapped to genename
  output_comb_type$genename_A = scer_annotation[orth_map_type_subset$target_genename, c('uniprot_id')]  #If data is mapped to uniprot id
  
  output_comb_type$LFC_A = outputA_all[output_comb_type$genename_A,output_cond]
  
  #Some s.cer orf names do not map to a protein name on the annotation file.  Set those values to 0.  
  no_protein_name_for_orf = which(is.na(output_comb_type$genename_A))
  output_comb_type$orth_type[no_protein_name_for_orf] = 'no_protein_name_for_orf'
  output_comb_type$LFC_A[no_protein_name_for_orf] = 0
  
  output_comb = rbind(output_comb, output_comb_type)
}


#Set expression values for all genes present in S. cer but not in other species
output_comb_type = data.frame(genename_A=setdiff(rownames(outputA_all), output_comb$genename_A), 
                              genename_B='NONE', 
                              LFC_B = 0, 
                              orth_type = 'Scer_only')

output_comb_type$LFC_A = outputA_all[output_comb_type$genename_A, output_cond]

scer_annotation_rlookup = scer_annotation
scer_annotation_rlookup = scer_annotation_rlookup[which(!(is.na(scer_annotation_rlookup$gene_name))),]
rownames(scer_annotation_rlookup) = scer_annotation_rlookup$gene_name

output_comb_type$sc_orf = scer_annotation_rlookup[output_comb_type$genename_A,"systematic_name"]
output_comb = rbind(output_comb, output_comb_type)



#plot scatter plots using Plotly 
lfc_scatter = ggplot(data=output_comb, aes(x=LFC_A, y=LFC_B, color=orth_type, text = paste(specA, ' name: ', genename_A, '\n', specB, ' name: ', genename_B, sep='')))  +   #key=genename_A
              geom_point(alpha=0.2)+
              labs(x=specA, y=specB, title=output_cond)

ggplotly(lfc_scatter)
#ggplotly(lfc_scatter,source = "select", tooltip = c("key"))


#summary(output_comb)
table(factor(output_comb$orth_type))
#table(factor(output_comb_type$orth_type))

