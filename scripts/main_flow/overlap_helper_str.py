# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 11:01:59 2022

@author: chenyibo
"""


import numpy as np
from itertools import permutations
import skbio
# from numba import jit


#%%

def generate_ovlp(reads,pen):
    edges = np.zeros(shape=(len(reads),len(reads)),dtype=np.float64)
    for i,j in permutations(range(len(reads)), 2):
        
        read1,read2 = reads[i],reads[j]
        res = skbio.alignment.local_pairwise_align_nucleotide(skbio.DNA(read1), skbio.DNA(read2), match_score=2, mismatch_score=-1)
        if res[1] > 10:
            align_block_len = res[2][0][1] - res[2][0][0] + 1
            mismatches = int((2*align_block_len - res[1])/4)

        score = -( align_block_len - pen * mismatches)

        edges[i, j] = score
    
    return edges

    



