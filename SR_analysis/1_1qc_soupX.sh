#!/bin/sh

out_dir=$1
in_dir=$2
sample_id=$3

echo "start soupX ${sample_id} in ${in_dir}"
Rscript --vanilla 1_1qc_soupX.R $out_dir $in_dir $sample_id
echo "finish soupX ${sample_id} in ${out_dir}"
