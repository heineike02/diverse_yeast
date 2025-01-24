import os
from evcouplings.utils import read_config_file
from evcouplings.utils.pipeline import execute

base_dir = os.path.normpath('/home/heineike/alphafold') #Ben's computer

print(base_dir)
config2 = read_config_file(base_dir + os.sep + os.path.normpath('examples/etc/coupling_analysis/ec_og1254_config.txt'))
print(config2)
outcfg = execute(**config2)