#Compile Ortholog Mapping between two species and make comparison for a given condition and two species

#Make sure to uncomment appropriate lines for working directories.

#base_dir = "/camp/home/heineib/working/Ben/"
#working_dir = paste(base_dir, 'diverse_strains/processed_data/', sep='' )
working_dir = "~/OneDrive - Charité - Universitätsmedizin Berlin/R_analysis/Proteomics/processed_data/"

#adding this for shits and giggles to test my understanding of Github workflow

#Compare LFC across species using ortholog mapping

output_cond = 'CN1_C2_2'

specA = 'Klac'
specB = 'Scer'


#Load Protein Annotation for S. cerevisiae which is used to map orf names onto S. cer data

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


#The New Uniprot id mapping doesn't have the double uniprot names listed and doesn't have NAs for Uniprot when a genename is present. 

#There's no way we could recover data for the genes that have NA for Uniprot but data linked to genename - I assume DIANN throws those out
#Not sure how diann deals with duplicate uniprot ids



#  Currently only set up to compare LFC (log fold change) since expression in proteomics experiments is relative so that is very hard to interpret.  

#Function to parse long name (used in uniprot annotation file and ortholog map) to get identifier for data.
# Example:  for sp|P05467|YKP1_KLULA, extracts P05467

orth_map_name_parse = function(name_split) {
  #name_out = strsplit(name_split[3], '_')[[1]][1] #in our example this would extract YKP1
  name_out = name_split[2]
  return(name_out)
}


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


output_comb = spec_compare(output_cond = output_cond, specA = specA, specB= specB, annotation_df= scer_annotation, working_dir=working_dir)



vis_option = 'sc_gal'

#Possible vis options: 
#orth_type:  Ortholog type 
#sc_gal: 4 top S. cerevisiae gal genes (defined below)
#pka_inh and pka_rep:  Pka repression genes from Heineike et al 2020. 

#  Gal Pathway 
# Main genes upregulated in S.cer from Dalal et al 2016
#P04385, GAL1_YEAST, YBR020W
#P04397, GAL10_YEAST, YBR019C 
#P08431, GAL7_YEAST, YBR018C 
#P13181, GAL2_YEAST, YLR081W 

sc_gal_list = c('P04385', 'P04397', 'P08431', 'P13181')

output_comb$sc_gal = FALSE
output_comb$sc_gal[which(output_comb$genename_B %in% intersect(output_comb$genename_B, sc_gal_list))] = TRUE


#  PKA genes
if (vis_option %in% c('pka_act', 'pka_rep')){
  sc_pka_inh = read.csv(file=paste(base_dir, 'external_data/pka_inhibition/20181017_deseq_SC_AS_WT_nmpp1.csv',sep=''), header=TRUE, row.names=1 )  
  
  #PKA activation and inhibition as defined in Heineike et al 2021
  pka_act = sc_pka_inh[which(sc_pka_inh$log2FoldChange>2.0 & sc_pka_inh$padj< 0.001),]
  pka_rep = sc_pka_inh[which( (sc_pka_inh$log2FoldChange< (-2.0)) & (sc_pka_inh$padj< 0.001)),]
  output_comb$pka_act= FALSE
  output_comb$pka_act[which(output_comb$sc_orf %in% intersect(output_comb$sc_orf, rownames(pka_act)))] = TRUE
  output_comb$pka_rep= FALSE
  output_comb$pka_rep[which(output_comb$sc_orf %in% intersect(output_comb$sc_orf, rownames(pka_rep)))] = TRUE
  }


# Other visualization options:  
#  Ohnologs
#  Cytoplasmic translation



#plot scatter plots using Plotly 
lfc_scatter = ggplot(data=output_comb, aes(x=LFC_A, y=LFC_B, color=get(vis_option), text = paste(get(vis_option), '\n',specA, ' name: ', genename_A, '\n', specB, ' name: ', genename_B, '\n Scer name: ', sc_name, sep='')))  +   #key=genename_A
  geom_point(alpha=0.2)+
  labs(x=specA, y=specB, title=output_cond)
lfc_scatter$labels$colour = vis_option

ggplotly(lfc_scatter, tooltip=c('x','y','text'))
a#ggplotly(lfc_scatter,source = "select", tooltip = c("key"))


#summary(output_comb)
table(factor(output_comb$orth_type))
#table(factor(output_comb_type$orth_type))
