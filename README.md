# diverse_yeast
 scripts for BMH's contributions to the Ralser Lab Yeast Structural Evolution Project.  Includes codes for 
   - Assembling orthogroup sequences
   - Making sequence based alignments
   - Refining structural alignment to prepare for DN/DS calculations
   - DN/DS Calculations

Core functions are contained in diverse yeast tools.py

Also uses a few functions from https://github.com/heineike02/y1000plus_tools
  make_og_genes_lookup
  extract_protein_seqs

Requires following directories: 
base_dir = os.path.normpath('G:/My Drive/Crick_LMS/projects/diverse_yeasts/alphafold')

y1000plus_dir = os.path.normpath('C:/Users/heineib/Documents/GitHub/y1000plus_tools/data') + os.sep

genomes_dir:  Includes genomic information from model organisms. 



# Key Scripts: 

## 20211025_alphafold_selection.ipynb
select orthogroups and orthologs based on selected species and annotations from yeast pathway database.  Builds list of proteins of interest based on orthogroups from orthomcl clusters

Visualize species trees for various figures. 

Assemble orthogroups and annotate metadata for each orthogroup in base_dir/selected_proteins/og_metadata.json.  Build a single fasta for all peptides in base_dir/selected_proteins/all_proteins.fasta. 

Quantify proteins that were calculated. 

20220719_alphafold_PLDDT_scores.ipynb

Assesses PLDDT scores of calculated structure and analyses structures that could not be calculated. 


20221012_assemble_peptide_cds_fastas.ipynb

Makes selected_protein_ids.json to list proteins we consider for each species. 

builds list of coding sequences for all proteins we consider for each species. 

Loads tm_align clustering data and assembles and assesses subclusters. 

Filters sequences for dN/dS analysis based on cluster size and sequence length.  

20220921_map_proteomes_to_shen.ipynb
Uses Pairwise sequence similarity scores to map uniprot ids to shen ids in two groups - first proteomes that were part of the proteomics experiment, then other proteomes with uniprot ids.  




20230724_phenotype_selection.ipynb
20240220_alignment_example.ipynb
20221206_selection_calculations.ipynb
20230201_ipath_plots.ipynb
20240312_thiamine_pathway.ipynb



# Proteome Annotation

Downloaded 'uniprot-proteome_UP000002311.tab' on 20220120 from uniprot after adding additional name fields to include the systematic name (Gene names  (ordered locus ))
