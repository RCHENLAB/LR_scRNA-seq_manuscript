#!/bin/bash

##### This script was modified based on sicelore
##### Please find original code and more information at https://github.com/ucagenomix/sicelore

# need Java.1.8 in JAVA_HOME
export PATH=./jdk1.8.0_121/bin/:$PATH
export JAVA_HOME=./jdk1.8.0_121/bin/java
export PATH=$PATH:./minimap2_10142020/minimap2:./racon/build/bin/racon:./poaV2/poa:./samtools-1.11/samtools
export PATH=$PATH:./minimap2_10142020/:./racon/build/bin/:./poaV2/:/storage/chen/home/jw29/software/samtools-1.11/

main="."
### BAM and barcode list from short read data
scbam="${main}/$1/10x/possorted_genome_bam"
cellbc="${main}/$1/10x/filtered_feature_bc_matrix/barcodes.tsv"

### Raw sequences (fastq) from long read data
gunzip -c ${main}/$1/$1.fastq.gz > ${main}/$1/$1.fq
fastq="${main}/$1/$1.fq"

shell_dir="./sicelore/Jar"
java=`which $JAVA_HOME`
poa=`which poa`
racon=`which racon`
minimap2=`which minimap2`
samtools=`which samtools`

output_dir="${main}/$1/sicelore"
tmp_dir="${output_dir}/tmp/"

junc_bed="./gencode.v25.mm10.junctions.bed"
fasta="./mm10.fasta.gz"
refFlat="./gencode.v25.mm10.refFlat.txt"

obj="${main}/$1/sicelore/10x_possorted_genome_bam.obj"
nanopore_read="$main/$1/sicelore/passed/$1FWD.fq"

if [ -z "$java" ] || [ -z "$poa" ] || [ -z "$samtools" ] || [ -z "$minimap2" ] || [ -z "$racon" ]
then
    echo -e "\nMissing path to required softwares:"
    echo -e "\tjava=$java"
    echo -e "\tpoa=$poa"
    echo -e "\tsamtools=$samtools"
    echo -e "\tminimap2=$minimap2"
    echo -e "\tracon=$racon"
    echo -e "\nPlease update your \$PATH and rerun.\n\n"
    exit
fi

# create output directory
mkdir $output_dir
mkdir $tmp_dir

## parse illumina bam file
$java -Xms30G -Xmx150G -XX:-UseGCOverheadLimit -XX:+UseConcMarkSweepGC  -jar  ${shell_dir}/IlluminaParser-1.0.jar -i ${scbam}.bam -o $obj -t ${cellbc} -b CB -g GN -u UB

# scan nanopore reads
$java -jar ${shell_dir}/NanoporeReadScanner-0.5.jar -i ${fastq} -o $output_dir
rm ${fastq}

# map reads to genome
$minimap2 -ax splice -uf --MD --sam-hit-only -t 6 --junc-bed $junc_bed $fasta $nanopore_read > $output_dir/minimap.sam
$samtools view -Sb $output_dir/minimap.sam -o $output_dir/minimap.unsorted.bam
$samtools sort $output_dir/minimap.unsorted.bam -o $output_dir/minimap.bam
$samtools index $output_dir/minimap.bam
rm $output_dir/minimap.sam $output_dir/minimap.unsorted.bam

# tag reads with gene name
$java -jar -Xmx4g ${shell_dir}/Sicelore-2.0.jar AddGeneNameTag I=$output_dir/minimap.bam O=$output_dir/GE.bam REFFLAT=$refFlat GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
$samtools index $output_dir/GE.bam

# tag reads with fastq sequence
$java -jar -Xmx250g -XX:-UseGCOverheadLimit -XX:+UseConcMarkSweepGC ${shell_dir}/Sicelore-2.0.jar AddBamReadSequenceTag I=$output_dir/GE.bam O=$output_dir/GEUS.bam FASTQ=$nanopore_read VALIDATION_STRINGENCY=SILENT
$samtools index $output_dir/GEUS.bam

# tag reads with cellBC/UMI barcodes
$java -jar -Xmx250g -XX:-UseGCOverheadLimit -XX:+UseConcMarkSweepGC ${shell_dir}/NanoporeBC_UMI_finder-1.0.jar -i $output_dir/GEUS.bam -o $output_dir/GEUS10xAttributes.bam -k $obj --ncpu 16 --maxUMIfalseMatchPercent 1 --maxBCfalseMatchPercent 5 --logFile $output_dir/cellBC_match.log
$samtools index $output_dir/GEUS10xAttributes.bam
$samtools index $output_dir/GEUS10xAttributes_umifound_.bam
