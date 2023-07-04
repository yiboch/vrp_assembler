#! /bin/bash

eval "$(conda shell.bash hook)"
conda activate bioinfo
rl=$1
for i in {1..10};do
    vn=$rl$i
    gfatools gfa2fa ./${rl}k/$i/v$vn.bp.hap1.p_ctg.gfa > ./${rl}k/$i/hap1.fasta
    gfatools gfa2fa ./${rl}k/$i/v$vn.bp.hap2.p_ctg.gfa > ./${rl}k/$i/hap2.fasta
    echo "Done converting sample $vn"
done
