#!/bin/sh

source anaconda3/bin/activate  scvi-env 

query=$1
celltype=$2
label=$3
query_name=$4

python 2label_transfer_general_sample.py $query $celltype $label $query_name 