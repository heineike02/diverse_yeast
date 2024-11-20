#!/bin/bash

python3 20241025_site_model_check_run.py


while read p; do
  echo "$p"
done < /home/heineike/alphafold_ln/selection_calculations/site_model/uncalculated_ogs.txt
