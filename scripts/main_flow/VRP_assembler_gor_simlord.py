# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 11:38:22 2021.

For solving VRP on Google OR-Tools, with Guided Local Search

@author: Lenovo-yibo
"""
#%%
# 生成 reads 和 pairwise alignment
import time
sst = time.time()
import numpy as np

from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp
import os
from importlib import reload
import overlap_helper_paf_simlord

print('[VRPasm]:Importing libraries complete, time cost:',round(time.time()-sst,3),'s')
st = time.time()
def print_solution(manager, routing, solution, num_reads):
    """Prints solution on console."""
    print(f'Objective: {solution.ObjectiveValue()}')
    routes = {}
    total_distance = 0
    SE = []
    for vehicle_id in range(2):
        tmp = []
        index = routing.Start(vehicle_id)
        plan_output = 'Route for vehicle {}:\n'.format(vehicle_id)
        route_distance = 0
        veh_range = []
        SEtmp = 0
        while not routing.IsEnd(index):
            plan_output += ' {} ->'.format(manager.IndexToNode(index))
            previous_index = index
            index = solution.Value(routing.NextVar(index))
            route_distance += routing.GetArcCostForVehicle(
                previous_index, index, vehicle_id)
            tmp.append((manager.IndexToNode(previous_index),manager.IndexToNode(index)))
            if manager.IndexToNode(index) == 1:
                veh_range = range(num_reads[0]+1)
            elif len(veh_range)!=0:
                prev_node = tmp[-1][0]
                this_node = tmp[-1][1]
                if (not this_node in veh_range) and (prev_node in veh_range):
                    SEtmp += 1
        routes[vehicle_id] = tmp
        plan_output += ' {}\n'.format(manager.IndexToNode(index))
        plan_output += 'Distance of the route: {}m\n'.format(route_distance)
        print(plan_output)
        total_distance += route_distance
        SE.append(SEtmp)
    print('Total Distance of all routes: {}m'.format(total_distance))
    print('Total switching:',max(SE))
    return routes

def check_overlap(overlap,i1,i2):
    print(overlap[i1,i2]+overlap[i2-1,i1+1],overlap[i1,i1+1]+overlap[i2-1,i2])

# dir_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/data/yeast/BFH/sim_reads/'
dir_path = "/mnt/d/projects/vrpassembler/data/yeast/BFH/sim_reads/"
paf_path = dir_path + 'aln.paf'

with open(dir_path+'nr.txt','r') as f:
    nr_str = f.readline()
    num_reads = list(map(int, nr_str.split(' ')))

#
# 开始优化
st = time.time()
reload(overlap_helper_paf_simlord)

pen1 = 5
pen2 = 0
print(f'[VRPasm]:Building overlap matrix with penalty {pen1} {pen2}', end='...')
overlap_, af_df, arigin_df = overlap_helper_paf_simlord.process_paf(paf_path, pen1, pen2)
overlap = overlap_.copy()

ovmin = overlap.min()
m1 = np.identity(overlap.shape[0]-1)
overlap[1:,1:] -= ovmin

m1 *= ovmin
overlap[1:,1:] += m1
et = time.time()
# if ovmin<0:
#     overlap[1:,1:] -= ovmin
    
#     m1 *= ovmin
#     overlap[1:,1:] += m1
# else:
#     ovmax = overlap.max()
#     tmp1 = overlap[1:,1:].copy()
#     tmp2 = np.where(tmp1 == 0, ovmax, tmp1)
#     m1 *= ovmax
#     tmp2 -= m1
#     overlap[1:,1:] = tmp2.copy()
print('DONE. time cost:',round(et-st,3),'s')

af_df['per_q'] = af_df['alignment_block_length'] / af_df['query_length']
af_df['per_t'] = af_df['alignment_block_length'] / af_df['target_length']

#%%
# 调用 ortools 
pd = False
pickups_deliveries = [(460,461),(535,536),(188,189)]

st = time.time()
print('[VRPasm]:Initialzing optimization', end='...')

dest_num = overlap.shape[0]

manager = pywrapcp.RoutingIndexManager(overlap.shape[0], 2, 0)

routing = pywrapcp.RoutingModel(manager)

def distance_callback(from_index, to_index):
    """Returns the distance between the two nodes."""
    # Convert from routing variable Index to distance matrix NodeIndex.
    from_node = manager.IndexToNode(from_index)
    to_node = manager.IndexToNode(to_index)
    return overlap[from_node,to_node]

transit_callback_index = routing.RegisterTransitCallback(distance_callback)

routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

if pd:
    for request in pickups_deliveries:
        p_id = manager.NodeToIndex(request[0])
        d_id = manager.NodeToIndex(request[1])
        # routing.AddPickupAndDelivery(p_id, d_id)
        # routing.AddPickupAndDelivery(d_id, p_id)
        routing.solver().Add(routing.VehicleVar(p_id) == routing.VehicleVar(d_id))

search_parameters = pywrapcp.DefaultRoutingSearchParameters()
search_parameters.first_solution_strategy = (
    routing_enums_pb2.FirstSolutionStrategy.LOCAL_CHEAPEST_INSERTION)
search_parameters.local_search_metaheuristic = (
    routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH)
search_parameters.time_limit.seconds = 300
search_parameters.log_search = True

print('DONE, time cost:', round(time.time()-st,3),'s')
st = time.time()
print('[VRPasm]:Start optimization')
print('='*80)
print('='*80)

solution = routing.SolveWithParameters(search_parameters)

print('[VRPasm]:Optimization DONE. time cost:',round(time.time()-st,3),'s')
st = time.time()

if solution:
    routes = print_solution(manager, routing, solution, num_reads)
    # switching_error = count_se(routes)
    
else:
    print('No solution found !')
    # breakpoint()

#%%
# routes = dict()
# base = 0
# for i in range(2):
#     num_r = num_reads[i]
#     routes[i] = [(0,base+1)]
#     routes[i].extend([(base+j, base+j+1) for j in range(1, num_r)])
#     routes[i].append((base+num_r, 0))
#     base = num_r
# 
# 提取路径，合成序列
print('='*80)
print('='*80)
import gc
rel = gc.collect()  # 回收多余变量的内存

from Bio import SeqIO
import time

reads_path = dir_path + 'real_reads.fasta'

st = time.time()
print('[VRPasm]:Indexing reads file', end='...')
index = SeqIO.index(reads_path, 'fasta')
et = time.time()
print('DONE. time cost', round(et-st,3),'s')

print('[VRPasm]:Reconstructing haplotypes', end='...')
haps = []
for veh, route in routes.items():
    
    reads_phase_file = dir_path + f'reads_{veh+1}.phase.fasta'
    with open(reads_phase_file,'w+') as ff:  # 存储分型好的reads
        for edge in route[:-1]:
            ff.write(f'>{edge[1]}')
            ff.write('\n')
            ff.write(str(index[f'{edge[1]}'].seq))
            ff.write('\n')

    hap_file = dir_path + f'hap{veh+1}.phase.fasta'
    i_hap = ''
    for edge in route[:-1]:
        if edge[0] == 0:
            i_hap += str(index[f'{edge[1]}'].seq)
        else:
            i_ov = af_df[(af_df['query_name']==edge[0]) & (af_df['target_name']==edge[1])]
            if i_ov.size == 0:
                i_ov = af_df[(af_df['query_name']==edge[1]) & (af_df['target_name']==edge[0])]
                if i_ov.size == 0:
                    i_hap += str(index[f'{edge[1]}'].seq)
                    continue
            if i_ov['query_name'].values[0]==edge[0]:
                backspace = i_ov['query_length'].values[0] - i_ov['query_start'].values[0]
                i_hap = i_hap[:-backspace] + str(index[f'{edge[1]}'].seq)[i_ov['target_start'].values[0]:]
            else:
                backspace = i_ov['target_length'].values[0] - i_ov['target_start'].values[0]
                i_hap = i_hap[:-backspace] + str(index[f'{edge[1]}'].seq)[i_ov['query_start'].values[0]:]

    haps.append(i_hap)
    with open(hap_file,'w+') as f:
        f.write(f'>hap{veh+1}')
        f.write('\n')
        f.write(i_hap)
        f.write('\n')

print('DONE. time cost',round(time.time()-et,3),'s\n')
print(f'length hap_1: {len(haps[0])}\nlength hap_2: {len(haps[1])}')

# %%
