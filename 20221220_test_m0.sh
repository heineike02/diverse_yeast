 #!/bin/bash

. /opt/conda/etc/profile.d/conda.sh
conda activate diverse_yeast_env

BASE_DIR=/home/heineike_wsl2/alphafold/
BASE_M0=${BASE_DIR}selection_calculations/m0/
 
OG_BASE=OG4150_REF_Scer_AF-P07256-F1-model_v2

CALC_DIR=${BASE_M0}${OG_BASE}/

cd $CALC_DIR

# Run codeML
/var/lib/paml/paml-4.10.6/bin/codeml ${CALC_DIR}m0.ctl
#ls ${CALC_DIR}