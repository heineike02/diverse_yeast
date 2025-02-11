import os
from evcouplings.utils import read_config_file
from evcouplings.utils.pipeline import execute

base_dir = os.path.normpath('/home/heineike/alphafold') #Ben's computer

#config = read_config_file(base_dir + os.sep + os.path.normpath('examples/etc/coupling_analysis/COX9_COX2/ec_config_COX9.txt'))
#config = read_config_file(base_dir + os.sep + os.path.normpath('examples/etc/coupling_analysis/COX9_COX2/ec_config_COX2.txt'))
config = read_config_file(base_dir + os.sep + os.path.normpath('examples/etc/coupling_analysis/COX9_COX2/ec_config_COX9_COX2_concat.txt'))
print(config)
outcfg = execute(**config)
