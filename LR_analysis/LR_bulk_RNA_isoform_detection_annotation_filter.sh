#!/bin/sh

#### Flair
cat ./fastq_pass/* > ./Ms_bulk/Ms_bulk.fastq.gz
pycoQC -f ./sequencing_summary.txt -o ./Ms_bulk/Ms_bulk_pycoQC.html

### flair align
flair align -g mm10.fasta -r ./Ms_bulk/Ms_bulk.fastq.gz -o ./Ms_bulk/flair_align --threads 36

### Correct misaligned splice sites using genome annotations
flair correct -q ./Ms_bulk/flair_align.bed -g mm10.fasta -f gencode.v25.mm10.annotation.gtf -o ./Ms_bulk/LR_corrected --threads 36

### Define high-confidence isoforms from corrected reads
flair collapse -q ./Ms_bulk/LR_corrected_all_corrected.bed -g mm10.fasta -f gencode.v25.mm10.annotation.gtf -t 36 -m /PATH/TO/minimap2 -sam /PATH/TO/samtools -o ./Ms_bulk/LR -r ./Ms_bulk/Ms_bulk.fastq.gz --stringent


#### SQANTI3
export PYTHONPATH=$PYTHONPATH:/PATH/TO/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/PATH/TO/cDNA_Cupcake/

main=./Ms_bulk
mkdir $main/SQANTI3

### QC
awk '$7!="."{print $0}' $main/LR.isoforms.gtf > $main/SQANTI3/Flair_LR_strand.gtf
python sqanti3_qc.py $main/SQANTI3/Flair_LR_strand.gtf gencode.v25.mm10.annotation.gtf mm10.fasta -d $main/SQANTI3 -o SQANTI3 --saturation --report both
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $main/SQANTI3/SQANTI3_corrected.gtf > $main/SQANTI3/LR_isoforms_corrected.gtf

### Filtering
python sqanti3_filter.py rules $main/SQANTI3/SQANTI3_classification.txt --gtf $main/SQANTI3/SQANTI3_corrected.gtf -d $main/SQANTI3 -o SQANTI3
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $main/SQANTI3/SQANTI3.filtered.gtf > $main/SQANTI3/LR_isoforms_filtered.gtf



