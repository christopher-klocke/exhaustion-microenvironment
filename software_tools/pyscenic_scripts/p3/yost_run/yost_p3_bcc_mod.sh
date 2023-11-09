echo "starting script"
# set this to directory that contains input file and will contain input and output files
run_directory=/Users/klockec/Documents/data/pyscenic_runs/yost_bcc_run1
code_directory=/Users/klockec/Documents/code/immune-tme-inference/software_tools/pyscenic_scripts
files_directory=/Users/klockec/Documents/data/pyscenic
run_filename=yost_bcc_pyscenic_input.loom

echo "running arboreto"
arboreto_with_multiprocessing.py $run_directory/$run_filename $files_directory/hs_hgnc_curated_tfs.txt \
    --method grnboost2 --output $run_directory/grnboost_out.tsv --num_workers 16 --seed 777

echo "running ctx"
pyscenic ctx $run_directory/grnboost_out.tsv \
    $files_directory/hg19-500bp-upstream-10species.mc9nr.feather \
    $files_directory/hg19-tss-centered-10kb-10species.mc9nr.feather \
    --annotations_fname $files_directory/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname $run_directory/$run_filename \
    --output $run_directory/ctx_output1.gmt \
    --mask_dropouts \
    --num_workers 16

echo "running aucell"
pyscenic aucell \
    $run_directory/$run_filename \
    $run_directory/ctx_output1.gmt \
    --output $run_directory/aucell_output1.loom \
    --num_workers 16

echo "script has finished."
