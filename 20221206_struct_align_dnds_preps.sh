#!/bin/bash
#source .venv_kits/bin/activate

. /opt/conda/etc/profile.d/conda.sh
conda activate diverse_yeast_env

cd /home/heineike_wsl2/github_s/diverse_yeast/

#python 20221206_struct_align_dnds_msas.py
python 20221206_struct_align_dnds_trees.py