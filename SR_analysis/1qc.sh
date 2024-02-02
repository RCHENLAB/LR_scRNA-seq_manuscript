#!/bin/sh

sample_id=""
in_dir=./${sample_id}/10x
main=./${sample_id}/SR
out_dir=${main}/SoupX
out_dir_rmDblt=${main}/rmDblt

mkdir -p $main
mkdir $out_dir

exec 1>>${main}/${sample_id}_log.txt
exec 2>>${main}/${sample_id}_log.txt

echo -e "${in_dir} ${sample_id} startTime:::::::\c" ; date
date_start=$(date +%s)
echo "## Pipeline Starts to Run##"

sh 1_1qc_soupX.sh $out_dir $in_dir $sample_id
sh 1_2qc_rmDoublet.sh $out_dir_rmDblt $out_dir $sample_id

echo -e "${in_dir} ${sample_id} endTime:::::::\c" ; date
date_end=$(date +%s)
echo "total time: $((date_end-date_start)) s"

