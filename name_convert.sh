#!/bin/bash
#activate environment for biokit 
. /home/heineike/.venv_biokit/bin/activate

BASE=/home/heineike/Crick_LMS/projects/diverse_yeasts/alphafold/selection_calculations/$1

biokit file_format_converter -i ${INPUT_BASE}${OG}_aln_clean.fasta -iff fasta -o ${OUTPUT_BASE}${OG}/${OG}_cds_aln.phy -off phylip

#deactivate environment for biokit 
deactivate