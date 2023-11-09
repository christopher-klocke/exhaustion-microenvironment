"""-run with 'sc2' conda env
-run with:
    python3 run2.py [run name (for log file)] [config file name]
-e.g.
    python3 run2.py tumor_tumor_overlap tumor_tumor_overlap.json
    python3 run2.py tumor_viral_overlap_li_wang tumor_viral_overlap_li_wang.json
-wrapper script for overlap analysis stage 1 and 2
"""
# setup logging
import sys
run_name = 'run_' + sys.argv[1] + '.log'
import logging
logging.basicConfig(filename=run_name, filemode='w', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# import statements
from overlap import overlap_analysis

# run analysis
with open(sys.argv[2], 'r') as input1:
    overlap_analysis(config_file=input1)
