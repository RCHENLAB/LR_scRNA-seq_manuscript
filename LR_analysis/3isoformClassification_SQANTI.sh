#!/bin/sh

##### This script was modified based on SQANTI3
##### Please find original code and more information at https://github.com/ConesaLab/SQANTI3

export PYTHONPATH=$PYTHONPATH:./cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:./cDNA_Cupcake/

main=./$1/Flair
mkdir $main/SQANTI3

### QC
awk '$7!="."{print $0}' $main/LR.isoforms.gtf > $main/SQANTI3/Flair_LR_strand.gtf
python sqanti3_qc.py $main/SQANTI3/Flair_LR_strand.gtf ./gencode.v25.mm10.annotation_nochr.gtf ./mm10.fasta -d $main/SQANTI3 -o SQANTI3 --saturation --report both
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $main/SQANTI3/SQANTI3_corrected.gtf > $main/SQANTI3/LR_isoforms_corrected.gtf

### Filtering
python sqanti3_filter.py rules $main/SQANTI3/SQANTI3_classification.txt --gtf $main/SQANTI3/SQANTI3_corrected.gtf -d $main/SQANTI3 -o SQANTI3
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $main/SQANTI3/SQANTI3.filtered.gtf > $main/SQANTI3/LR_isoforms_filtered.gtf

### Rescue
python sqanti3_rescue.py rules $main/SQANTI3/SQANTI3_RulesFilter_result_classification.txt --isoforms $main/SQANTI3/SQANTI3_corrected.fasta -k $main/SQANTI3/SQANTI3_classification.txt --gtf $main/SQANTI3/SQANTI3.filtered.gtf -f ./mm10.fasta -g ./gencode.v25.mm10.annotation_nochr.gtf -d $main/SQANTI3 -o SQANTI3 -j utilities/filter/filter_default.json
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $main/SQANTI3/SQANTI3_rescued.gtf > $main/SQANTI3/LR_isoforms_rescued.gtf
