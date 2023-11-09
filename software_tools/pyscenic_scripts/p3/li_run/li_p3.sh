echo "starting script"
# set this to directory that contains input file and will contain input and output files
run_directory=/Users/klockec/Documents/data/analysis_files/p3/pyscenic_runs/p3/li_run
run_filename=li_pyscenic_input.loom

echo "running arboreto"
python /Users/klockec/Documents/code/immune-tme-inference/software_tools/pyscenic_scripts/arboreto_run.py 20 $run_directory/grnboost_out.tsv \
    /Users/klockec/Documents/data/pyscenic/hs_hgnc_curated_tfs.txt \
    $run_directory/$run_filename

echo "running ctx"
python /Users/klockec/Documents/code/immune-tme-inference/software_tools/pyscenic_scripts/custom_start.py ctx \
    --mode custom_multiprocessing --num_workers 20 \
    --expression_mtx_fname $run_directory/$run_filename \
    -o $run_directory/ctx_output1.gmt \
    --annotations_fname /Users/klockec/Documents/data/pyscenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    $run_directory/grnboost_out.tsv \
    /Users/klockec/Documents/data/pyscenic/hg19-500bp-upstream-10species.mc9nr.feather \
    /Users/klockec/Documents/data/pyscenic/hg19-tss-centered-10kb-10species.mc9nr.feather

echo "running aucell"
python /Users/klockec/Documents/code/immune-tme-inference/software_tools/pyscenic_scripts/custom_start.py aucell \
    --num_workers 20 \
    -o $run_directory/aucell_output1.loom \
    $run_directory/$run_filename \
    $run_directory/ctx_output1.gmt

echo "script has finished."
