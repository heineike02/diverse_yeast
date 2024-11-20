#!/bin/bash

#Needed to mount google drive before running this
#sudo mount -t drvfs G: /mnt/g  

MITO_PDB_BASE=/home/heineike/alphafold_ln/examples/etc/mitochondrial_sequences/pdbs

for MITO_CORE_GENE in cox1 cox2 cox3 cob atp6 atp8 atp9
do
	PDB_DIR=${MITO_PDB_BASE}/pdbs/
	PDB_LIST=${MITO_PDB_BASE}/pdb_list_${MITO_CORE_GENE}.txt
	/home/heineike/usalign/USalign -dir $PDB_DIR $PDB_LIST -suffix .pdb -mm 4 > ${MITO_PDB_BASE}/us_align_${MITO_CORE_GENE}.fasta
done
