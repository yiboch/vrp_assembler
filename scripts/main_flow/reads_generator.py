# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 17:27:09 2022

@author: chenyibo
"""
import random
import numpy as np


def generate_reads(len_avg_rds, sig_rds_len, len_dna, overlap, sig_ovlp, snp_rate, 
                   num_ploid, len_repeat, num_repeats, error_rate, path, error_on_snp):
    
    bases = ['A','T','G','C']
    choices = random.choices
    sample = random.sample
    ori_dna_seq = choices(bases,k=len_dna)
    all_seqs = [ori_dna_seq.copy() for _ in range(num_ploid)]

    repeat_seq = choices(bases,k=len_repeat)
    ins_sites = choices(range(len_dna),k=num_repeats)
    ins_sites.sort()
    ins_after = []
    delta = -len_repeat
    for i in ins_sites:
        delta += len_repeat
        ii = i + delta
        ins_after.append(ii)
        for j in range(num_ploid):
            all_seqs[j] = np.append(np.append(all_seqs[j][:ii], repeat_seq), all_seqs[j][ii:])
    
    
    
    n_snps = round(snp_rate*len_dna)
    sites = sample(range(len_dna), k=n_snps)
    
    for site in sites:
        tmp_site = []
        if num_ploid>4:
            tmp_site = choices(bases,k=num_ploid)
        else:
            tmp_site = sample(bases,k=num_ploid)
        for i in range(num_ploid):
            all_seqs[i][site] = tmp_site[i]

    all_reads = [[] for i in range(num_ploid)]

    
    reads_info = []
    ii = -1
    
    for i in range(num_ploid):
        
        ioverlap = round(max(5,random.gauss(mu=overlap, sigma=sig_ovlp)))
        ilen = ioverlap + max(1,round(random.gauss(mu=len_avg_rds-ioverlap, sigma=sig_rds_len)))
        reach = ilen
        while True:
            if reach < len_dna + num_repeats*len_repeat:
            
                iread = all_seqs[i][reach-ilen : reach].copy()
                ii+=1
                for j in range(len(iread)):

                    if (reach-ilen+j in sites) and (not error_on_snp):
                        continue
                    
                    reads_info.append([ii,j])
                        
                all_reads[i].append(iread)
                ioverlap = round(max(5,random.gauss(mu=overlap, sigma=sig_ovlp)))
                while ioverlap>=ilen:
                    ioverlap = round(max(5,random.gauss(mu=overlap, sigma=sig_ovlp)))
                ilen = ioverlap + max(1, round(random.gauss(mu=len_avg_rds-ioverlap, sigma=sig_rds_len)))
                
                reach += (ilen - ioverlap)
                
                
            else:
                iread = all_seqs[i][max(reach-ilen,ioverlap + 1) : ].copy()
                ii+=1

                for j in range(len(iread)):
                    
                    if (reach-ilen+j in sites) and (not error_on_snp):
                        continue
                    
                    reads_info.append([ii,j])
                        
                all_reads[i].append(iread)
                break
        
        
    reads = []

    for i in range(num_ploid):
        reads.extend(all_reads[i])
        
    deno_n_bases = len(reads_info)
    n_error = round(error_rate*deno_n_bases)
    
    err_sites = sample(range(deno_n_bases), k = n_error)
    
    for ies in err_sites:
        rec = reads_info[ies]
        iread = reads[rec[0]]
        
        ii = rec[1]
        
        erid = choices([1,2,3,4],k=1)[0]

        # point mutation
        if erid==1:#
            nbases = bases.copy()
            nbases.remove(iread[ii])
            iread[ii] = choices(nbases,k=1)[0]
        elif erid==2:
            iread[ii] = ''
        elif erid==3:
            nbases = bases.copy()
            nbases.remove(iread[ii])
            iread[ii] += choices(nbases,k=1)[0]
        else:
            iread[ii] *= 2
    
    read_file = path + rf'/{num_ploid}_{len_dna+num_repeats*len_repeat}_reads.fasta'
    with open(read_file, mode='w+') as f:
        i = 1
        for rd in reads:
            f.write(f'>read_{i}\n')
            f.write(''.join(rd))
            f.write('\n')
            i += 1
    
    return reads, read_file
            
if __name__=='__main__':
    
    len_avg_rds = 15
    len_dna = 50
    overlap = 10
    avg_snp_dis = 15
    num_ploid = 2
    len_repeat = 0
    num_repeats = 0
    
    data_path = 'D:/projects/SNP/data/'
    
    reads = generate_reads(len_avg_rds, len_dna, overlap, avg_snp_dis, 
                   num_ploid, len_repeat, num_repeats, data_path)
    