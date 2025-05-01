##Functions for use in mapping proteome data between yeast orthologs

library(ggplot2)
library(GGally)
library(plotly)

base_dir = "/camp/home/heineib/working/Ben/"
working_dir = paste(base_dir, 'diverse_strains/processed_data/', sep='' )
#working_dir = "~/OneDrive - Charité - Universitätsmedizin Berlin/R_analysis/Proteomics/processed_data/"

cat(paste("base_dir: ", base_dir, '\nworking_dir: ', working_dir, '\nEnsure base_dir and working_dir are set correctly for your setup in diverse_yeast_tools.R\n'))


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



#Function to parse long name (used in uniprot annotation file and ortholog map) to get identifier for data.
# Example:  for sp|P05467|YKP1_KLULA, extracts P05467

orth_map_name_parse = function(name_split) {
  #name_out = strsplit(name_split[3], '_')[[1]][1] #in our example this would extract YKP1
  name_out = name_split[2]
  return(name_out)
}






#Combines data from two species.  Currently the default specB is S.cer and it is not set up for other species.
spec_compare = function(output_cond, specA, specB= 'Scer', annotation_df, working_dir) {
  #Function inputs:  
  #  output_cond:  Output condition that is being compared - must be in the columns of the output file
  #  specA:  First species in the comparison - 
  #  specB:  Second species in the comparison default is S. cerevisiae BY4741
  #  annotation_df: Dataframe used for annotation.  Default is scer_annotation.  
  #
  #  This function requires the file <specA>_<specB>.csv to be saved in <working_dir>#processed_data/ortholog_maps/
  #
  #  Returns:  output_comb, a dataframe with lines for each gene and the following fields: 
  #   genename_<X> where X is A or B.  If specB is S. cerevisiae the names are the Uniprot Ids.
  #   <output>_<X> where output is the output_type ('LFC' or 'exp')
  #   sc_orf:  systematic name for S. cerevisiae.  Required for ortholog mapping when specB is S. cerevisiae
  #   sc_name:  Standard name for S. cerevisiae genes.  More easily interpreted.   
  
  #annotation_df = scer_annotation
  
  outputA_all = read.csv(file=paste(working_dir , specA, '/LFC_data_',specA,'.csv', sep=''), header=TRUE, row.names=1) #output_list[[output_type]][[specA]] 
  
  #Mapping S.cerevisae to lab strain, BY4741
  if (specB=='Scer') {
    outputB_all = read.csv(file=paste(working_dir, 'Scer-BY4741KI/LFC_data_Scer-BY4741KI.csv', sep=''), header=TRUE, row.names=1) 
  } else {
    outputB_all = read.csv(file=paste(working_dir, specB, '/LFC_data_',specB,'.csv', sep=''), header=TRUE, row.names=1) 
  }
  
  
  
  #Load Ortholog Mapping
  orth_map = read.csv(file=paste(working_dir, 'ortholog_maps/', specA, '_', specB, '.csv', sep=''), header=TRUE, row.names=1 )  
  
  orth_map_source_genename = strsplit(orth_map$source_genename, '[|]')
  
  orth_map_source_new_names = sapply(orth_map_source_genename, orth_map_name_parse)
  
  orth_map$source_genename_short = orth_map_source_new_names
  
  #Check all output genes are included in ortholog map
  missing_from_orthmap = setdiff(rownames(outputA_all), orth_map$source_genename_short)
  if (length(missing_from_orthmap>0)) {
    warning(paste(length(missing_from_orthmap)), ' genes missing from ortholog map for ', specA, sep='')
  }
  
  output_comb = data.frame()
  
  #No orthologs in specB
  #Sets default LFC to 0
  for (orth_type in c('no_eggnog_orthologs','no_target_orthologs' )) {
    orth_map_type_subset = orth_map[which(orth_map$orth_type==orth_type),]
    output_comb_type = data.frame(genename_A = orth_map_type_subset$source_genename_short, 
                                  LFC_A = outputA_all[orth_map_type_subset$source_genename_short, output_cond],
                                  genename_B = 'NONE', 
                                  sc_orf = 'NONE', 
                                  sc_name = 'NONE',
                                  LFC_B = 0, 
                                  orth_type = orth_type
    )
    
    output_comb = rbind(output_comb, output_comb_type)
  }
  
  #orthologs present in specB
  for (orth_type in c('one2one', 'no_target_orthologs', 'one2many', 'many2one', 'many2many')) {
    orth_map_type_subset = orth_map[which(orth_map$orth_type==orth_type),]
    output_comb_type = data.frame(genename_A = orth_map_type_subset$source_genename_short, 
                                  LFC_A = outputA_all[orth_map_type_subset$source_genename_short, output_cond],
                                  orth_type = orth_type, 
                                  sc_orf = orth_map_type_subset$target_genename, 
                                  sc_name = 'NONE', #Default sc_name
                                  LFC_B= NA #Default LFC
    )
    
    output_comb_type$sc_name = scer_annotation[orth_map_type_subset$target_genename, c('Entry.name')]  
    output_comb_type$genename_B = scer_annotation[orth_map_type_subset$target_genename, c('uniprot_id')]  
    
    output_comb_type$LFC_B = outputB_all[output_comb_type$genename_B,output_cond]
    
    
    
    #Genes from specB that do not map to the annotation file.  Set those values to 0.  
    no_protein_name_for_orf = which(is.na(output_comb_type$genename_B))
    output_comb_type$orth_type[no_protein_name_for_orf] = 'no_protein_name_for_orf'
    output_comb_type$LFC_B[no_protein_name_for_orf] = 0
    output_comb = rbind(output_comb, output_comb_type)
  }
  
  
  #Set expression values for all genes present in SpecB but not in SpecA
  output_comb_type = data.frame(genename_B=setdiff(rownames(outputB_all), output_comb$genename_B), 
                                genename_A='NONE', 
                                LFC_A = 0, 
                                orth_type = paste(specB,'_only', sep='')
  )
  
  output_comb_type$LFC_B = outputB_all[output_comb_type$genename_B, output_cond]
  
  scer_annotation_rlookup = scer_annotation
  scer_annotation_rlookup = scer_annotation_rlookup[which(!(is.na(scer_annotation_rlookup$gene_name))),]
  rownames(scer_annotation_rlookup) = scer_annotation_rlookup$gene_name
  
  output_comb_type$sc_orf = scer_annotation_rlookup[output_comb_type$genename_B,"systematic_name"]
  output_comb_type$sc_name = scer_annotation_rlookup[output_comb_type$genename_B,"Entry.name"]
  output_comb = rbind(output_comb, output_comb_type)
  
  
  return(output_comb)
}


load_scer_annotation = function(){
  #Loads scerevisiae annotation file

  scer_annotation = read.delim(file=paste(working_dir, 'uniprot-proteome_UP000002311.tab', sep=''),sep='\t', header=TRUE)

  #Switched to uniprot annotation file because it was missing some assignments from uniprot gene name to uniprot ID
  # 
  #Previously used annotations_proteins.tsv which was from the dataframe annotations_proteins from 210318_annotation_tables.Rdata
  #write.table(annotations_proteins, paste(working_dir,'annotations_proteins.tsv', sep = ''))
  
  #renames columns to match column names from older annotation file
  scer_annotation = rename(scer_annotation, c('uniprot_id'='Entry', 'systematic_name'='Gene.names...ordered.locus..'))
  
  scer_all_data = read.csv(file=paste(working_dir, 'Scer-BY4741KI/LFC_data_Scer-BY4741KI.csv', sep=''), header=TRUE, row.names=1) 
  print("List of Uniprot IDs present in the data that are not in the annotation file:")
  print(setdiff(rownames(scer_all_data), scer_annotation$uniprot_id))
  rownames(scer_annotation) = scer_annotation$systematic_name
  
  #Assess items in the annotation file that have duplicate uniprot ids and are NA for uniprot ID, but have a genename assigned.
  uniprot_split = strsplit(scer_annotation$uniprot_id, ',')
  uniprot_lengths = lapply(uniprot_split,length)
  scer_annotation_uniprot_doubles = scer_annotation[which(uniprot_lengths==2),]
  scer_annotation_uniprot_na = scer_annotation[which(is.na(scer_annotation$uniprot_id)),]
  scer_annotation_uniprot_na_genename_present = scer_annotation_uniprot_na[which(!is.na(scer_annotation_uniprot_na$gene_name)),]
  #write.table(scer_annotation_uniprot_doubles, paste(working_dir, 'scer_duplicate_uniprot_ids.tsv', sep=''))
  #write.table(scer_annotation_uniprot_na_genename_present, paste(working_dir, 'scer_uniprot_na_genename_present.tsv', sep=''))
  
  
  #The Uniprot_id mapping doesn't have the double uniprot names listed and doesn't have NAs for Uniprot when a genename is present. 
  
  #There's no way we could recover data for the genes that have NA for Uniprot but data linked to genename - I assume DIANN throws those out
  #Not sure how diann deals with duplicate uniprot ids
  
  return(scer_annotation)
}



data_table_from_scer_list = function(output_type, conds, other_spec_list){
  #output data frame:  
  
  #
  
  # Scer_systematic_name, 
  # Scer_uniprot_id, 
  # Scer_genename, 
  # Scer_<cond>,   # Output value of S.cer for he given condition
  # other_spec, 
  # ortholog_name, 
  # orth_type,
  # <cond>
  
  
  sc_map = data.frame()
  
  for (spec in other_spec_list) {
    #spec = other_spec_list[1]
    
    #Load species data
    output_spec = read.csv(file=paste(working_dir, spec, '/', output_type, '_data_',spec,'.csv', sep=''), header=TRUE, row.names=1) 
    
    #Load Ortholog Mapping
    orth_map = read.csv(file=paste(working_dir, 'ortholog_maps/', spec, '_',source_spec ,'.csv', sep=''), header=TRUE, row.names=1 )  
    
    orth_map_source_genename = strsplit(orth_map$source_genename, '[|]')
    
    orth_map_source_new_names = sapply(orth_map_source_genename, orth_map_name_parse)
    
    orth_map$source_genename_short = orth_map_source_new_names
    
    #Check all output genes are included in ortholog map
    missing_from_orthmap = setdiff(rownames(output_spec), orth_map$source_genename_short)
    if (length(missing_from_orthmap>0)) {
      warning(paste(length(missing_from_orthmap)), ' genes missing from ortholog map for ', spec, sep='')
      warning(paste(missing_from_orthmap, collapse=' '))
    }
    
    #Map orthologs from S.cerevisiae systematic names onto the species. 
    sc_map_spec_raw = orth_map[which(orth_map$target_genename %in% sc_gene_list), c('target_genename', 'source_genename_short', 'orth_type')]
    
    #Adds S.cerevisiae orthologs of e.g. one2many genes for which only one paralog was provided. 
    
    sc_map_spec = orth_map[which(orth_map$source_genename_short %in% sc_map_spec_raw$source_genename_short), c('target_genename', 'source_genename_short', 'orth_type')]
    
    if (length(sc_map_spec$target_genename) > length(sc_map_spec_raw$target_genename)) {
      warning(paste('The following S. cerevisae genes were not in the original list and were added due to orthology with genes from ', spec, sep =''))
      warning(paste(setdiff(sc_map_spec$target_genename, sc_map_spec_raw$target_genename), collapse = ' '))
    }
    
    sc_map_spec$other_spec = spec
    
    #Map data from original species
    sc_map_spec = cbind(sc_map_spec, output_spec[sc_map_spec$source_genename_short,conds])
    
    sc_map = rbind(sc_map, sc_map_spec)
  }
  
  #warnings()
  
  #Map on S_cer annotations columns 'uniprot_id' and 'Entry.name'
  scer_annotation = load_scer_annotation()
  sc_map = merge(sc_map, scer_annotation[,c('uniprot_id', 'systematic_name','Entry.name')], by.x='target_genename', by.y='systematic_name', all.x=TRUE, all.y=FALSE)
  
  #Map on S_cer data
  output_scer = read.csv(file=paste(working_dir, 'Scer-BY4741KI/', output_type, '_data_Scer-BY4741KI.csv', sep=''), header=TRUE, row.names=1)
  colnames(output_scer) =  paste('Scer_', colnames(output_scer))
  sc_map = cbind(sc_map, output_scer[sc_map$uniprot_id,])
  
  return(sc_map)
  
  
}

