#!/bin/sh

##### This script was modified based on Flair
##### Please find original code and more information at https://flair.readthedocs.io/en/latest/index.html

mkdir ./Flair

### Extract SR junctions and filter out junctions with few supporting junction reads [3]
#junctions_from_sam -s $1/10x/possorted_genome_bam.bam -n $1/Flair/junctions_from_SR
#awk '{ if ($5 > 3) { print } }' $1/Flair/junctions_from_SR_junctions.bed | awk '{gsub(/^chr/,""); print}'> $1/Flair/junctions_from_SR.filtered.bed

### Prepare LR BED12 
bedtools bamtobed -bed12 -i ./sicelore/demultiplexed.bam > ./Flair/LR.bed

### Correct misaligned splice sites using genome annotations and short-read splice junctions
flair correct -q ./Flair/LR.bed -g ./mm10.fasta -f ./gencode.v25.mm10.annotation_nochr.gtf -o ./Flair/LR_corrected
#flair correct -q ./Flair/LR.bed -g ./mm10.fasta -f ./gencode.v25.mm10.annotation_nochr.gtf -j ./Flair/junctions_from_SR.filtered.bed -o ./Flair/LR_corrected

### Define high-confidence isoforms from corrected reads
cat ./sample1/sicelore/passed/sample1FWD.fq ./sample2/sicelore/passed/sample2FWD.fq ./sample3/sicelore/passed/sample3FWD.fq ./sample4/sicelore/passed/sample4*FWD.fq > ./Flair/reads.fq
flair collapse -q ./Flair/LR_corrected_all_corrected.bed -g ./mm10.fasta -f ./gencode.v25.mm10.annotation_nochr.gtf -t 36 -m ./minimap2_10142020/minimap2 -sam ./samtools-1.16.1/samtools -o ./Flair/LR -r ./Flair/reads.fq --stringent

### Extract LR read name with cellBC
less ./sicelore/combined.sam | awk -F '\t' '{ split($0,a,"BC:Z:"); split(a[2],b,"\t"); if($1!~"@") print $1"\t"b[1]}' > ./Flair/LR_ID_BC
sed -i '1iread_id\tbarcode' ./Flair/LR_ID_BC
python3 ./barcode.py ./Flair/LR_ID_BC ./SR_barcode_celltype ./Flair/LR_ID_celltype

### Quantify FLAIR isoform usage across samples using minimap2
awk '{ if($3=="AC") print $1}' ./Flair/LR_ID_celltype > ./Flair/AC
awk '{ if($3=="BC") print $1}' ./Flair/LR_ID_celltype > ./Flair/BC
awk '{ if($3=="Cone") print $1}' ./Flair/LR_ID_celltype > ./Flair/Cone
awk '{ if($3=="MG") print $1}' ./Flair/LR_ID_celltype > ./Flair/MG
awk '{ if($3=="Rod") print $1}' ./Flair/LR_ID_celltype > ./Flair/Rod
awk '{ if($3=="RGC") print $1}' ./Flair/LR_ID_celltype > ./Flair/RGC
awk '{ if($3=="BC1A") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/BC1A
awk '{ if($3=="BC1B") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/BC1B
awk '{ if($3=="BC2") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/BC2
awk '{ if($3=="BC3A") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/BC3A
awk '{ if($3=="BC3B") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/BC3B
awk '{ if($3=="BC4") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/BC4
awk '{ if($3=="BC5A") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/BC5A
awk '{ if($3=="BC5B") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/BC5B
awk '{ if($3=="BC5C") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/BC5C
awk '{ if($3=="BC5D") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/BC5D
awk '{ if($3=="BC6") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/BC6
awk '{ if($3=="BC7") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/BC7
awk '{ if($3=="BC8") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/BC8
awk '{ if($3=="BC9") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/BC9
awk '{ if($3=="RBC") print $1}' ./Flair/LR_ID_celltype_BC > ./Flair/RBC
awk '{ if($3=="Glycinergic") print $1}' ./Flair/LR_ID_celltype_AC > ./Flair/Glycinergic
awk '{ if($3=="GABAergic") print $1}' ./Flair/LR_ID_celltype_AC > ./Flair/GABAergic
awk '{ if($3=="nGnG") print $1}' ./Flair/LR_ID_celltype_AC > ./Flair/nGnG

for x in AC BC Cone MG Rod RGC BC1A BC1B BC2 BC3A BC3B BC4 BC5A BC5B BC5C BC5D BC6 BC7 BC8 BC9 RBC Glycinergic GABAergic nGnG
do
seqtk subseq ./Flair/reads.fq ./Flair/$x > ./Flair/$x.fq
rm ./Flair/$x
done

flair quantify -r ./Flair/reads_manifest.tsv -i ./Flair/LR.isoforms.fa -t 36 -m ./minimap2_10142020/minimap2 -sam ./samtools-1.16.1/samtools -o ./Flair/LR.flair.quantify --generate_map --tpm --stringent


### Plots
plot_isoform_usage ./Flair/LR.isoforms.bed ./Flair/LR.flair.quantify ENSMUSG00000000617

