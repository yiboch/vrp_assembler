#! /bin/bash

eval "$(conda shell.bash hook)"
conda activate minimap2

hap_filename=$1
ref_base='../data/references/MHC_37'
echo "File name is $hap_filename"
base_filename=${hap_filename%.*}

minimap2 -t 12 -a --secondary=no $ref_base.mmi $hap_filename | samtools view -bS | samtools sort > $base_filename.aln.sort.bam
samtools index $base_filename.aln.sort.bam
samtools mpileup -g -f $ref_base.fa $base_filename.aln.sort.bam | bcftools filter -i 'DP=1' -O b | bcftools call -v -m > $base_filename.phase.vcf


echo "dealing phase file, now process indexing"
bgzip $base_filename.phase.vcf
tabix -0 -p vcf $base_filename.phase.vcf.gz


