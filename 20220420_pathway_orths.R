source('~/github/diverse_yeast/diverse_yeast_tools.R')

# Load the package required to read JSON files.
library("rjson")
library("ape")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("treeio")
BiocManager::install(version='3.14')
BiocManager::install("ggtree")
#BiocManager::install("ComplexHeatmap")

library("treeio")
library("ggtree")

library(ggplot2)
library(tidyr)
library(tibble)
#library(hrbrthemes)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
#library(patchwork)

#Goal:  To show the presence or absence of pathway orthologs for each pathway.  

#Make a named list for each pathway that contains the following data: 

#sc_genes
#orthologs - a matrix with species and the number of orthologs for each of the sc_genes
#orthologs_eggnog - a matrix with species and the number of orthologs using my eggnog ortholog mapping
#orthologs_expression
#orthologs_eggnog_expression


# Give the input file name to the function.
#Note:  Needed to replace NaN with "NULL" in the json file to conform to json standards not enforced in python.  Better to fix the way the file is saved. 
og_metadata <- fromJSON(file = paste(base_dir, 'diverse_strains/alphafold/og_metadata_nanrep.json', sep=''), unexpected.escape = 'skip')

#Build a dataframe that maps yeast genes to orthogroups
sc_gene_og_map = data.frame()
for (og_name in names(og_metadata)) {
  og = og_metadata[[og_name]]
  sc_gene_og_map = rbind(sc_gene_og_map, data.frame(og$sc_genes, og_name))
}

row.names(sc_gene_og_map) = sc_gene_og_map[,1]

#Load pathways table
pathways = read.csv(paste(base_dir, 'diverse_strains/alphafold/pathway_list.tsv', sep=''), sep = '\t')

#Load species table
species_table = read.csv(paste(base_dir, 'diverse_strains/alphafold/species_selection.csv', sep=''))
species_ind_map = species_table[which(species_table['Load']=='Y'),c('Species.name', 'spec_og_id', 'Time_tree_name')]


#Import time calibrated tree and order species table by tree

time_tree = read.newick(paste(base_dir, "diverse_strains/332_2408OGs_time-calibrated_phylogeny_species-names_updated.newick", sep=""))

time_tree_subset <- drop.tip(time_tree, setdiff( time_tree$tip.label, species_ind_map$Time_tree_name))

tree_plot = ggplot(time_tree_subset, aes(x, y)) + geom_tree() + theme_tree()  #geom_tiplab()

spec_order = get_taxa_name(tree_plot)

tree_plot



spec_ind_fun = function(og_gene) {
  as.numeric(strsplit(og_gene, '_')[[1]][1])
}


orth_check_fun = function(spec_ind_check, spec_inds) {
  sum(spec_inds == spec_ind_check)
}


pathway_orth_map_list = list()
pathway_min_orth_present = list()
pathway_pct_orth_missing = list()
pathway_n_genes = list()


#for each pathway
for (pathway_genes_raw in pathways[['sc_genes']]) {
  
  pathway_split = strsplit(pathway_genes_raw, "'")[[1]]
  
  pathway_genes = pathway_split[seq(from=2,to=length(pathway_split)-1,by=2)]
  
  
  #Cycle through genes in the pathway and make a matrix that tells whether they are present in each of the species. 
  
  #Input: species subset
  
  orth_presence_df = data.frame(row.names = species_ind_map[['Time_tree_name']])
  
  for (sc_gene in pathway_genes) {
    #sc_gene = pathway_genes[1]
    
    #Find its orthogroup
    og_name = sc_gene_og_map[sc_gene,'og_name']
    
    og_genes = og_metadata[[og_name]]$og_genes
    
    
    #See which species are present in that orthogroup.  
    
    spec_inds = sapply(og_genes, spec_ind_fun)
    
    orth_presence_df[sc_gene] = sapply(species_ind_map[['spec_og_id']], orth_check_fun, spec_inds)
  }
  
  
  
  
  #row.names(orth_presence_df) = orth_presence_df$row_names
  
  #orth_presence_df = orth_presence_df[,2:length(orth_presence_df)]
  
  if (length(pathway_genes)==1) {
    orth_presence_pct = c(sum(orth_presence_df>0)/length(rownames(orth_presence_df)))
  } else {
    orth_presence_pct = colSums(orth_presence_df>0)/length(rownames(orth_presence_df))
    
  }
  
  
  
  pathway_orth_map_list = append(pathway_orth_map_list, list(orth_presence_df))
  
  #calculate min % of genes present for any ortholog in the pathway
  pathway_min_orth_present = append(pathway_min_orth_present,min(orth_presence_pct))
  
  
  #calculate % genes in the pathway less than 100%
  pathway_pct_orth_missing = append(pathway_pct_orth_missing, sum(orth_presence_pct<1)/length(orth_presence_pct))
  
  #add length of pathway
  pathway_n_genes = append(pathway_n_genes, length(pathway_genes))
  
  
}


pathways$orth_map = pathway_orth_map_list

pathways$min_orth_present = unlist(pathway_min_orth_present)
pathways$pct_orth_missing = unlist(pathway_pct_orth_missing)
pathways$n_genes = unlist(pathway_n_genes)


p <- ggplot(pathways , aes(min_orth_present, pct_orth_missing, label = pathway_name))+
  geom_point(alpha=0.2) 
#geom_text(data=subset(pathways, min_orth_present <0.50),
#           aes(min_orth_present, pct_orth_missing, label = pathway_name), nudge_y = 0.02, hjust = 0)
#geom_text()

p

a = pathways[,c('pathway_name', 'min_orth_present', 'pct_orth_missing', 'n_genes')]

a[which((a['n_genes']>0) & (a['pct_orth_missing']>0.7)  &  (a['min_orth_present'] >0)),]


#full pathways: 
#L-lysine biosynthesis IV 
#pentose phosphate pathway (non-oxidative branch)

#pathway missing a few genes
#de novo biosynthesis of purine nucleotides

#pathways missing many genes
#'superpathway of NAD biosynthesis'
#'riboflavin, FMN and FAD biosynthesis'
#'aerobic respiration, electron transport chain'

ex_pathway = 'aerobic respiration, electron transport chain' #'riboflavin, FMN and FAD biosynthesis'# 'superpathway of NAD biosynthesis'#'de novo biosynthesis of purine nucleotides'# 'pentose phosphate pathway (non-oxidative branch)'#'L-lysine biosynthesis IV'  #'glutathione biosynthesis' #'biotin biosynthesis' # 'glutathione-glutaredoxin redox reactions'

pathways[which(pathways['pathway_name']==ex_pathway), 'sc_genes_name']

ex_orth_map = pathways[[which(pathways['pathway_name']==ex_pathway), 'orth_map']]

sc_genes_name_raw = pathways[which(pathways['pathway_name']==ex_pathway), 'sc_genes_name']

sc_genes_name_split = strsplit(sc_genes_name_raw, "'")[[1]]

sc_genes_name = sc_genes_name_split[seq(from=2,to=length(sc_genes_name_split)-1,by=2)]


colnames(ex_orth_map) = sc_genes_name

#heatmap(as.matrix(ex_orth_map), keep.dendro = FALSE, scale='none')


col_fun = colorRamp2(c(0, 1, max(as.matrix(ex_orth_map))), c("grey", "yellow", "red"))
lgd_at =seq(0, max(as.matrix(ex_orth_map)), by = 1)
hmap = Heatmap(as.matrix(ex_orth_map), col = col_fun,
               row_order = spec_order,  
               row_names_side = 'left', 
               row_names_max_width = max_text_width(rownames(as.matrix(ex_orth_map))), 
               heatmap_legend_param = list(
                 at = lgd_at,
                 title = "N genes", 
                 legend_gp = gpar(fill = col_fun(lgd_at)), 
                 title_position = "leftcenter-rot")
)




hmap



#Tried to combine heatmap and dendrogram - will just do it separately
#library(cowplot)
#library(grid)

#a = grid.grabExpr(hmap)

#plot_grid(tree_plot, grid.grabExpr(draw(hmap)), nrow=1, ncol =2)

#tree_plot + hmap



a = subset( pathways[,c("pathway_name","min_orth_present", "pct_orth_missing")], min_orth_present <0.50)

b = subset( pathways[,c("pathway_name","min_orth_present", "pct_orth_missing")], pct_orth_missing >0.80)

a# Heatmap 
ex_orth_map %>%
  
  # Data wrangling
  as_tibble() %>%
  rowid_to_column(var="X") %>%
  gather(key="Y", value="Z", -1) %>%
  
  ## Change Y to numeric
  mutate(Y=as.numeric(gsub("V","",Y)))   -> a  #%>%
#mutate(X=as.numeric(gsub("V","",X))) %>%

# Viz
#ggplot(aes(X, Y, fill= Z)) + 
#geom_tile() 
#theme_ipsum() +
#theme(legend.position="none")



#+ scale_size_continuous(range=c(0,30)) + #scale_area()+
#geom_point(aes(colour = oc$percent_women)) + 
coord_equal() +
  scale_colour_gradient(high = "red")+
  ylim(700, 1700) +
  xlim(700, 1700) +
  geom_abline(slope=1) +
  labs(title = "Income Disparity by Occupation and Gender") +
  ylab("Women's Weekly Earnings in $") +
  xlab("Men's Weekly Earnings in $")




