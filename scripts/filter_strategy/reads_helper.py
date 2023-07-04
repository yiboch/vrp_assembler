# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 17:27:09 2022

@author: chenyibo
"""
# import os
from Bio import SeqIO

def get_reads(fa_dir):
    
    ori_fasta_file = fa_dir + 'sim_reads_tmp.fa'
    fasta_file = fa_dir + 'sim_reads.fasta'

    reads = []

    for iread in SeqIO.parse(ori_fasta_file,'fasta'):
        reads.append(str(iread.seq))
    
    with open(fasta_file,'w+') as ff:
        i = 1
        for read in reads:
            ff.write(f'>{i}')
            ff.write('\n')
            ff.write(read)
            ff.write('\n')
            i+=1

    