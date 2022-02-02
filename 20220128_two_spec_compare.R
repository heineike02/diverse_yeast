source('~/github/diverse_yeast/diverse_yeast_tools.R')

#Compile Ortholog Mapping between two species and make comparison for a given condition and two species

#Make sure to uncomment appropriate lines for working directories.

#adding this for shits and giggles to test my understanding of Github workflow

#Compare LFC across species using ortholog mapping

output_cond = 'CN1_C2_2'

specA = 'Calb'
specB = 'Scer'


#Load Protein Annotation for S. cerevisiae which is used to map orf names onto S. cer data
scer_annotation = load_scer_annotation()

#  Currently only set up to compare LFC (log fold change) since expression in proteomics experiments is relative so that is very hard to interpret.  
output_comb = spec_compare(output_cond = output_cond, specA = specA, specB= specB, annotation_df= scer_annotation, working_dir=working_dir)



vis_option = 'gluc_rep'

#Possible vis options: 
#orth_type:  Ortholog type 
#sc_gal: 4 top S. cerevisiae gal genes (defined below)
#gluc_rep: glucose repression genes from Kayikci_Nielsen et al ????
#pka_inh and pka_rep:  Pka repression genes from Heineike et al 2020. 

#  Gal Pathway 
# Main genes upregulated in S.cer from Dalal et al 2016
#P04385, GAL1_YEAST, YBR020W
#P04397, GAL10_YEAST, YBR019C 
#P08431, GAL7_YEAST, YBR018C 
#P13181, GAL2_YEAST, YLR081W 

sc_gal_list = c('P04385', 'P04397', 'P08431', 'P13181')


#Gal genes from Dalal et al 2016

if (vis_option %in% c('gal_dalal_sc','gal_dalal_ca')) {
  gal_dalal_sc = read.csv(file=paste(base_dir, 'diverse_strains/processed_data/external_data/gal_dalal/sc_gal.csv', sep=''))
  output_comb$gal_dalal_sc= FALSE
  output_comb$gal_dalal_sc[which(output_comb$sc_orf %in% intersect(output_comb$sc_orf, gal_dalal_sc$Systematic.Name))] = TRUE
}





output_comb$sc_gal = FALSE
output_comb$sc_gal[which(output_comb$genename_B %in% intersect(output_comb$genename_B, sc_gal_list))] = TRUE


#  PKA genes
if (vis_option %in% c('pka_act', 'pka_rep')){
  sc_pka_inh = read.csv(file=paste(base_dir, 'diverse_strains/processed_data/external_data/pka_inhibition/20181017_deseq_SC_AS_WT_nmpp1.csv',sep=''), header=TRUE, row.names=1 )  
  
  #PKA activation and inhibition as defined in Heineike et al 2021
  pka_act = sc_pka_inh[which(sc_pka_inh$log2FoldChange>2.0 & sc_pka_inh$padj< 0.001),]
  pka_rep = sc_pka_inh[which( (sc_pka_inh$log2FoldChange< (-2.0)) & (sc_pka_inh$padj< 0.001)),]
  output_comb$pka_act= FALSE
  output_comb$pka_act[which(output_comb$sc_orf %in% intersect(output_comb$sc_orf, rownames(pka_act)))] = TRUE
  output_comb$pka_rep= FALSE
  output_comb$pka_rep[which(output_comb$sc_orf %in% intersect(output_comb$sc_orf, rownames(pka_rep)))] = TRUE
  }


# Gluc Repression
if (vis_option == 'gluc_rep') {
  gluc_rep = read.csv(file=paste(base_dir, 'diverse_strains/processed_data/external_data/gluc_rep/Proteins_C_metabolite_repression_Kayikci_Nielsen_final.csv', sep=''))
  output_comb$gluc_rep= FALSE
  output_comb$gluc_rep[which(output_comb$genename_B %in% intersect(output_comb$genename_B, gluc_rep$uniprot_id))] = TRUE
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
#ggplotly(lfc_scatter,source = "select", tooltip = c("key"))


#summary(output_comb)
table(factor(output_comb$orth_type))
#table(factor(output_comb_type$orth_type))
