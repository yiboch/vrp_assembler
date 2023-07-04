#! /bin/bash

eval "$(conda shell.bash hook)"
conda activate minimap2

hap_filename=$1
reads_filename=$2
echo "File name is $hap_filename"
base_filename=${hap_filename%.*}

echo "Base file name is $base_filename"

minimap2 -t 12 --secondary=no $hap_filename $reads_filename > $base_filename.racon.gfa.paf
racon -t 12 -w 2000 -e 0.002 $reads_filename $base_filename.racon.gfa.paf $hap_filename > $base_filename.racon.fasta
