#!/bin/bash

. /opt/conda/etc/profile.d/conda.sh
conda activate diverse_yeast_env

BASE_DIR=/home/heineikeb/alphafold/
#/home/heineike_wsl2/alphafold/
BASE_M1=${BASE_DIR}selection_calculations/m1/


#for OG_BASE in OG4150_REF_Scer_AF-P07256-F1-model_v2 OG2603_REF_Scer_AF-P50076-F1-model_v2 OG2845_REF_Scer_AF-P43577-F1-model_v2 OG3677_REF_Scer_AF-P47125-F1-model_v2 OG1299_REF_Scer_AF-P00549-F1-model_v2

#OG4*.tm.fasta.clipkit.renamed
#${BASE_DIR}msas/structural/tm_align/trees/
#for A in /home/heineikeb/alphafold/msas/structural/tm_align/trees/OG4*.tm.fasta.clipkit.treefile.renamed; 
#do 
#    echo $A
#done
#for TREE_FILE in /home/heineikeb/alphafold/msas/structural/tm_align/trees/OG*.tm.fasta.clipkit.treefile.renamed
#for OG_BASE in OG1299_REF_Scer_AF-P00549-F1-model_v2 OG1390_REF_Scer_AF-P47164-F1-model_v2 

#OG4552_REF_Scer_AF-P37299-F1-model_v2
#OG5441_REF_Scer_AF-P00724-F1-model_v2 OG1122_REF_Scer_AF-P13711-F1-model_v2 OG4150_REF_Scer_AF-P07256-F1-model_v2 OG4316_REF_Scer_AF-P00424-F1-model_v2 OG4346_REF_Scer_AF-P10174-F1-model_v2 OG4360_REF_Scer_AF-P08525-F1-model_v2 OG1377_REF_Scer_AF-P53090-F1-model_v2 OG2112_REF_Scer_AF-P33421-F1-model_v2 OG1631_REF_Scer_AF-P22217-F1-model_v2 OG4744_REF_Scer_AF-P22289-F1-model_v2
#OG4552_REF_Scer_AF-P37299-F1-model_v2 OG4555_REF_Scer_AF-P32799-F1-model_v2
#OG4755_kluyveromyces_marxianus__OG4755__51_759 OG5006_kluyveromyces_marxianus__OG5006__51_4218
#OG1275_REF_Scer_AF-P47052-F1-model_v2 OG1287_REF_Scer_AF-P37298-F1-model_v2 OG2006_REF_Scer_AF-P07257-F1-model_v2 OG2248_REF_Scer_AF-P07143-F1-model_v2 OG2704_REF_Scer_AF-P21801-F1-model_v2 OG2714_REF_Scer_AF-P00427-F1-model_v2 OG3208_REF_Scer_AF-P00128-F1-model_v2 OG3505_REF_Scer_AF-P04037-F1-model_v2 OG3575_REF_Scer_AF-P08067-F1-model_v2 OG4118_REF_Scer_AF-P07255-F1-model_v2 OG4352_REF_Scer_AF-P00127-F1-model_v2 OG4751_REF_Scer_AF-P04039-F1-model_v2 OG5490_REF_Scer_AF-P32340-F1-model_v2

for OG_BASE in OG2197_REF_Scer_AF-P00942-F1-model_v2 OG1145_REF_Scer_AF-P00358-F1-model_v2 OG1218_REF_Scer_AF-P08417-F1-model_v2 
do 
    #DIR_OG_BASE=$(echo $TREE_FILE | cut -d '.' -f 1)
    #OG_BASE=$(echo $DIR_OG_BASE | cut -d '/' -f 9)
    #echo $TREE_FILE
    #echo $DIR_OG_BASE
    echo $OG_BASE
    #Make directory
    CALC_DIR=${BASE_M1}${OG_BASE}/
    echo $CALC_DIR    
    mkdir $CALC_DIR
    
    #Change into directory and call codeML
    cd $CALC_DIR

    ## Copy in tree
    cp ${BASE_DIR}msas/structural/tm_align/trees/${OG_BASE}.tm.fasta.clipkit.treefile.renamed ${CALC_DIR}tree.treefile
    
    # Copy in phy
    cp ${BASE_DIR}msas/structural/tm_align/cds_trim_strict/${OG_BASE}.tm.fasta.clipkit.cds.renamed.codeML.phy ${CALC_DIR}aln.phy
    
    # Copy in control file
    cp ${BASE_M1}m1.ctl ${CALC_DIR}m1.ctl
    
    # Run codeML
    /var/lib/paml/paml-4.10.6/bin/codeml ${CALC_DIR}m1.ctl


    # remove data
    # for now skipping it
    # paml_results_out = base_dir + os.sep + os.path.normpath('selection_calculations/yn00/' + og_base + '.yn00.out')
    # shutil.copy(calc_dir + og_base + '.yn00.out', paml_results_out)


    #python /home/heineike_wsl2/github_s/diverse_yeast/20221206_struct_align_yn00.py

done
