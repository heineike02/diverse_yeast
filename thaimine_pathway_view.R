# #Install pathview
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("pathview")

library(pathview)
data(gse16873.d)
data(demo.paths)
i = 1
pway = demo.paths$sel.paths[i]

pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = pway, species = "hsa", out.suffix = "gse16873", kegg.native = T)

#pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110", species = "hsa", out.suffix = "gse16873")

#Check organism id
#data(korg)
#head(korg)
#
#S cerevisiae is 	
#ktax.id = T00005
#tax.id = 4932
#kegg.code = sce

thiamine_pway = "00730"
spec = "sce"

thia_data= read.csv('G:/My Drive/Crick_LMS/projects/diverse_yeasts/alphafold/examples/thiamine/test_data.txt', sep = "\t", row.names = 1)

pv.out <- pathview(gene.data=thia_data, pathway.id = thiamine_pway, species = spec, gene.idtype="kegg", out.suffix = "thiamine", kegg.native = T)



