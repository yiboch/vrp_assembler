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
    # random.seed(seed)
    choices = random.choices
    sample = random.sample
    # rand01 = np.random.random

    # 生成序列，生成SNP位点
    ori_dna_seq = choices(bases,k=len_dna)
    all_seqs = [ori_dna_seq.copy() for _ in range(num_ploid)]
    
    
    # 插入重复序列
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
            # breakpoint()
            all_seqs[j] = np.append(np.append(all_seqs[j][:ii], repeat_seq), all_seqs[j][ii:])
    
    
    
    n_snps = round(snp_rate*len_dna)
    sites = sample(range(len_dna), k=n_snps)
    
    # 改成均匀采样的部分
    for site in sites:
        tmp_site = []
        if num_ploid>4:
            tmp_site = choices(bases,k=num_ploid)
        else:
            tmp_site = sample(bases,k=num_ploid)
        for i in range(num_ploid):
            all_seqs[i][site] = tmp_site[i]

    # snp_site = random.randint(0, avg_snp_dis)
    # snps = []
    # while snp_site < len_dna:
    #     tmp_site = []
    #     if num_ploid>4:
    #         tmp_site = choices(bases,size=num_ploid)
    #     else:
    #         tmp_site = choices(bases,size=num_ploid,replace=False)
    #     for i in range(num_ploid):
    #         all_seqs[i][snp_site] = tmp_site[i]
    #     snps.append([snp_site, *tmp_site])
    #     sites.append(snp_site)
    #     snp_site += random.randint(int(avg_snp_dis/2), avg_snp_dis)  #round(random.gauss(mu=avg_snp_dis, sigma=5))
    
    
    # 生成reads
    all_reads = [[] for i in range(num_ploid)]

    
    reads_info = []
    ii = -1
    
    # ilen = len_avg_rds
    
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
    # error_str = []

    for i in range(num_ploid):
        reads.extend(all_reads[i])
        
    # 在分母的 error 个数
    deno_n_bases = len(reads_info)
    n_error = round(error_rate*deno_n_bases)
    
    err_sites = sample(range(deno_n_bases), k = n_error)
    
    for ies in err_sites:
        rec = reads_info[ies]
        iread = reads[rec[0]]
        
        ii = rec[1]
        
        # ierror = [rec[0],iread.copy()]
        
        erid = choices([1,2,3,4],k=1)[0] # 1234分别代表替换、删除、插入、复制

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
        
        # ierror.append(iread.copy())
        # error_str.append(ierror)
    
    # 保存文件
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
    