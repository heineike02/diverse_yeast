#Plot scatter plots of LFC for KL vs SC, ZR vs SC, and KL vs SC. 

base_dir = "/camp/home/heineib/working/Ben/"
working_dir = paste(base_dir, 'diverse_strains/processed_data/', sep='' )

#Visualize 
#  Ohnologs
#  Cytoplasmic translation
#  PKA genes
sc_pka_inh = read.csv(file=paste(base_dir, 'external_data/pka_inhibition/20181017_deseq_SC_AS_WT_nmpp1.csv',sep=''), header=TRUE, row.names=1 )  

pka_act = sc_pka_inh[which(sc_pka_inh$log2FoldChange>2.0 & sc_pka_inh$padj< 0.001),]
pka_rep = sc_pka_inh[which( (sc_pka_inh$log2FoldChange< (-2.0)) & (sc_pka_inh$padj< 0.001)),]
#  Gal Pathway 
# Main genes upregulated in S.cer from Dalal et al 2016
#P04385, GAL1_YEAST, YBR020W
#P04397, GAL10_YEAST, YBR019C 
#P08431, GAL7_YEAST, YBR018C 
#P13181, GAL2_YEAST, YLR081W 

sc_gal_list = c('P04385', 'P04397', 'P08431', 'P13181')



#Compare different LFC comparisons across species using ortholog mapping

#Start with 'CN1_C2_2'
#Tried 'CN1_2_

output_cond = 'CN1_C2_2'
specA = 'Scer'
specB = 'Zrou'
output_type = 'LFC'

outputA_all = read.csv(file=paste(working_dir , specA, '/', output_type, '_data_',specA,'.csv', sep=''), header=TRUE, row.names=1) #output_list[[output_type]][[specA]] 
outputB_all = read.csv(file=paste(working_dir, specB, '/', output_type, '_data_',specB,'.csv', sep=''), header=TRUE, row.names=1) #output_list[[output_type]][[specB]] 

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
print(setdiff(rownames(outputA_all), scer_annotation$uniprot_id))

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
#write.table(scer_annotation_uniprot_doubles, paste(working_dir, 'scer_duplicate_uniprot_ids.tsv', sep=''))
#write.table(scer_annotation_uniprot_na_genename_present, paste(working_dir, 'scer_uniprot_na_genename_present.tsv', sep=''))


#The New Uniprot id mapping doesn't have the double uniprot names listed and doesn't have NAs for Uniprot when a genename is present. 

#There's no way we could recover data for the genes that have NA for Uniprot but data linked to genename - I assume DIANN throws those out
#Not sure how diann deals with duplicate uniprot ids



output_comb = data.frame()
#colnames(orth_map_comb) = c('genename_B', 'LFC_B', 'genename_A', 'sc_orf','sc_name', 'LFC_A', 'orth_type')

for (orth_type in c('no_eggnog_orthologs','no_target_orthologs' )) {
  orth_map_type_subset = orth_map[which(orth_map$orth_type==orth_type),]
  output_comb_type = data.frame(genename_B = orth_map_type_subset$source_genename_short, 
                                LFC_B = outputB_all[orth_map_type_subset$source_genename_short, output_cond],
                                genename_A = 'NONE', 
                                sc_orf = 'NONE', 
                                sc_name = 'NONE',
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
                                sc_name = 'NONE',
                                LFC_A = NA
  )
  
  output_comb_type$sc_name = scer_annotation[orth_map_type_subset$target_genename, c('Entry.name')]  #If data is mapped to genename
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
output_comb_type$sc_name = scer_annotation_rlookup[output_comb_type$genename_A,"Entry.name"]
output_comb = rbind(output_comb, output_comb_type)


output_comb$pka_act= FALSE
output_comb$pka_act[which(output_comb$sc_orf %in% intersect(output_comb$sc_orf, rownames(pka_act)))] = TRUE
output_comb$pka_rep= FALSE
output_comb$pka_rep[which(output_comb$sc_orf %in% intersect(output_comb$sc_orf, rownames(pka_rep)))] = TRUE

output_comb$sc_gal = FALSE
output_comb$sc_gal[which(output_comb$genename_A %in% intersect(output_comb$genename_A, sc_gal_list))] = TRUE


#color=orth_type
#plot scatter plots using Plotly 
lfc_scatter = ggplot(data=output_comb, aes(x=LFC_A, y=LFC_B, color=sc_gal, text = paste(specA, ' name: ', genename_A, '\n', specB, ' name: ', genename_B, '\n Scer name: ', sc_name, sep='')))  +   #key=genename_A
  geom_point(alpha=0.2)+
  labs(x=specA, y=specB, title=output_cond)

ggplotly(lfc_scatter)
#ggplotly(lfc_scatter,source = "select", tooltip = c("key"))


#summary(output_comb)
table(factor(output_comb$orth_type))
#table(factor(output_comb_type$orth_type))

