#!/bin/bash

. /opt/conda/etc/profile.d/conda.sh
conda activate diverse_yeast_env
#conda activate diverse_yeast_env_ete3

python /home/heineikeb/github/diverse_yeast/20240916_struct_align_dnds_msas_core_surface.py
