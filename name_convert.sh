#!/bin/bash
#activate environment for biokit 
. /home/heineike/.venv_biokit/bin/activate
#deactivate with the command
#
#deactivate 
#
#when complete

YSE_BASE=/home/heineike/Crick_LMS/projects/diverse_yeasts/alphafold/
INPUT_BASE=${YSE_BASE}msas/ogs_cds_threaded_short/
OUTPUT_BASE=${YSE_BASE}selection_calculations/

for OG in OG1299 OG1390 
    do 
    OUTPUT_DIR=${OUTPUT_BASE}${OG}
    if [ ! -d "$OUTPUT_DIR" ]; then
        mkdir $OUTPUT_DIR
    fi
    biokit file_format_converter -i ${INPUT_BASE}${OG}_aln_clean.fasta -iff fasta -o ${OUTPUT_BASE}${OG}/${OG}_cds_aln.phy -off phylip
done 

#deactivate environment for biokit 
deactivate