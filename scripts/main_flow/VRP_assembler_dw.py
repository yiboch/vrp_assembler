#%%
# import libraries
import sys
import os
project_dir = '/mnt/d/projects/vrpassembler/'
# dw_dir = '/mnt/d/dw/D-Wave-VRP-master'
sys.path.append(os.path.join(project_dir, 'dw'))
sys.path.append(project_dir)
sys.path.append(project_dir+'codes_bak/')
from importlib import reload
import vrp_solvers
import vrp_problem
from input import *
import pandas as pd
import numpy as np
from overlap_helper_str import generate_ovlp
from reads_generator_str_bak import generate_reads
import Levenshtein

def get_qubo_energy(sample, qubo):
    qubo_energy = 0
    for k,v in qubo.dict.items():
        x1, x2 = k[0], k[1]
        qubo_energy += sample[x1] * sample[x2] * v
    return qubo_energy

#
# setup parameters
nn = 1

len_avg_rds = int(20*nn)
sig_rds_len = 3

len_dna = int(50*nn)

overlap = int(15*nn)
sig_ovlp = 0

avg_snp_dis = min(overlap-1,1000)
num_ploid = 3
len_repeat = 0
num_repeats = 0

error_rate = 0.0

eos = False

pen = 20

max_error_n = round(error_rate*len_avg_rds*num_ploid*(1+(len_dna+num_repeats*len_repeat-len_avg_rds)/overlap))

data_path = '/mnt/d/projects/vrpassembler/data/'
#
#  generate reads
reads, all_seqs, error_info, sites = generate_reads(len_avg_rds, sig_rds_len, len_dna, overlap, sig_ovlp, avg_snp_dis, 
                                                        num_ploid, len_repeat, num_repeats, 
                                                        error_rate, data_path, max_error_n, eos)
    
    
nbases = 0
for i in reads:
    nbases += len(i)
error_rate_actual = error_info[0]/nbases
# reads_nb = np.array(reads)
overlap = generate_ovlp(reads,pen)

#
# run annealer

# overlap = pd.read_csv(project_dir+'/data/overlap.csv',index_col=[0])
# overlap = overlap.to_numpy()

dest_num = overlap.shape[0]

costs = np.zeros((dest_num+1,dest_num+1))
costs[1:,1:] = overlap.copy()

capacities = [1 for _ in range(num_ploid)] # n cars

sources = [0]

dests = list(np.arange(dest_num)+1)

# reads 长度，这里简单先设为 20-mers
reads_len = [0]*(dest_num + 1)

#%%
reload(vrp_problem)
reload(vrp_solvers)

only_one_const = 25
order_const = 1.
problem = vrp_problem.VRPProblem(sources, costs, capacities, dests, reads_len)

solver = vrp_solvers.FullQuboSolver(problem)
qubo = problem.get_full_qubo(only_one_const, order_const)
maxqubo = max(max(qubo.dict.values()),-min(qubo.dict.values()))
print('maximum qubo strength', maxqubo)
#%%
import datetime
rr = list(range(5,100,5))
rr.insert(0, 0.5)
#%%
tmp_res = []
# for ant in rr:
ant = 0.5
now = str(datetime.datetime.now().time())
solution, sample = solver.solve(qubo, solver_type = 'hybrid_1', label=f'VRPasm_{now}', 
                        num_reads=200, auto_scale=True, annealing_time=ant, 
                        chain_strength=maxqubo+5)
#
# qubo_energy = get_qubo_energy(sample, qubo)
# Checking if solution is correct.
ttc = np.nan
if solution == None or solution.check() == False:
    print("Solution : ", solution.solution) 
    ttc = solution.total_cost()
    print("Total cost : ", ttc)
    print("Solver hasn't find solution or the solution is infeasible.\n")
else:
    print("Solution : ", solution.solution) 
    ttc = solution.total_cost()
    print("Total cost : ", ttc)
tmp_res.append(ttc)

# print("qubo energy : ", qubo_energy)
# print("\n")
print(tmp_res)
#%%
routes = solution.solution
hd = Levenshtein.hamming
# 提取路径
res = []
seqs = []
for i in routes:
    # breakpoint()
    if i == []:
        continue
    
    edges = {}
    for ii in range(len(i)-1):
        edges[i[ii]] = i[ii+1]
    
    res.append(edges)
    
    cur_seq = 0
    next_seq = edges[cur_seq]
    iseq = ''
    
    while next_seq != 0:
        read2 = reads[next_seq - 1]
        if iseq == '':
            iseq = read2
            print(iseq, ' ', '_', len(read2), next_seq )
            # print('_', len(read2), next_seq - 1)
        else:
            read1 = reads[cur_seq - 1]
            iovlp = overlap[cur_seq-1,next_seq-1]
            min_flip = min(len(read1),len(read2))
            for k in range(5, min_flip)[::-1]:
                l = hd(read1[-k:], read2[:k])
                # l = np.count_nonzero(read1[-min(k,len(read1)):], read2[:min(k,len(read2))])
                if - (k - pen * l) == iovlp:
                    break

            print(' '*(len(iseq)-k)+read2, ' ', k, len(read2), next_seq)
            # print(k, len(read2), next_seq-1)
            iseq += read2[k:]
        cur_seq = next_seq
        next_seq = edges[cur_seq]
    
    seqs.append(iseq)
    print(iseq)
    print('-'*len(iseq))

#%%
# 比对得分
reseq1 = seqs[0]

oriseqs = [''.join(all_seqs[i]) for i in range(num_ploid)]

alignments = np.array([Levenshtein.distance(oriseqs[i],reseq1) for i in range(num_ploid)])

match_ind = alignments.argmin()
alignment = alignments[match_ind]
matchseq = oriseqs[match_ind]