#!/bin/sh

out_dir=$1
in_dir=$2
sample_id=$3

echo "start DoubletFinder ${sample_id} in ${in_dir}"
Rscript --vanilla 1_2qc_rmDoublet.R $in_dir $sample_id $out_dir
echo "finish DoubletFinder ${sample_id} in ${out_dir}"