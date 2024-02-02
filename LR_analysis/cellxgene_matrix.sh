#!/bin/sh

less $1.sam | awk -F '\t' '{ split($0,x,"GE:Z:"); split(x[2],y,"\t"); split($0,a,"BC:Z:"); split(a[2],b,"\t"); if($1!~"@") print y[1]"\t"b[1]}' > $1.sicelore.mtx
less $1.sicelore.mtx | awk '{ print $1"\t"$2}' | sort | uniq -c | awk '{ print $2"\t"$3"\t"$1}' > $1.sicelore.mtx

Rscript Ms_cell_gene.R $1
