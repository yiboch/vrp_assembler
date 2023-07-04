#! /bin/bash

eval "$(conda shell.bash hook)"
conda activate simlord

avg_len=$1
min_len=$2
vn=$3
dir="../data/v$vn"
echo "Dir name is $dir"
if [ ! -d $dir ]; then
    mkdir $dir
    echo "Creating $dir"
fi

simlord -rr H1.fa -ln 0.02 0 $avg_len -mr $min_len -c 50 -mp 1 -pi 0.001 -ps 0.001 -pd 0.001 $dir/r1 
simlord -rr H2.fa -ln 0.02 0 $avg_len -mr $min_len -c 50 -mp 1 -pi 0.001 -ps 0.001 -pd 0.001 $dir/r2 

samtools view -b $dir/r1.sam | samtools sort > $dir/r1.sort.bam
samtools index $dir/r1.sort.bam

samtools view -b $dir/r2.sam | samtools sort > $dir/r2.sort.bam
samtools index $dir/r2.sort.bam  

python3 ./filter_strategy/pipeline.py -vn $vn

minimap2 -c -k 28 --for-only -D -t 12 -x ava-pb $dir/sim_reads.fasta $dir/sim_reads.fasta > $dir/aln.paf