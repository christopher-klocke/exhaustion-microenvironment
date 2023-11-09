"""-run with 'sc2' conda env
-run with:
    python3 run.py [run name (for log file)] [config file name]
-e.g.
    python3 run.py dev_run1 config_dev.json
-wraps extremes_differential step and bh correction into one wrapper script
-NOTE: input .h5ad file must be formatted so that (uniformly named) cell type annotations are present
  in adata.obs['celltype']
"""
# setup logging
import sys
run_name = 'run_' + sys.argv[1] + '.log'
import logging
logging.basicConfig(filename=run_name, filemode='w', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# import statements
from compare import analyze

# run analysis
with open(sys.argv[2], 'r') as input1:
    analyze(config_file=input1)
