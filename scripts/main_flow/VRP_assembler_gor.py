# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 11:38:22 2021.

For solving VRP on Gurobi, with DFJ formulation implemented

@author: Lenovo-yibo
"""
#%%
# import libraries
from reads_generator import generate_reads
import overlap_helper_paf
import numpy as np
from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp

from importlib import reload
import subprocess
import time
#
def print_solution(data, manager, routing, solution):
    """Prints solution on console."""
    print(f'Objective: {solution.ObjectiveValue()}\n')
    routes = {}
    # total_distance = 0
    for vehicle_id in range(data['num_vehicles']):
        tmp = []
        index = routing.Start(vehicle_id)
        plan_output = 'Route for vehicle {}:\n'.format(vehicle_id)
        # route_distance = 0
        while not routing.IsEnd(index):
            plan_output += ' {} ->'.format(manager.IndexToNode(index))
            previous_index = index
            index = solution.Value(routing.NextVar(index))
            # route_distance += routing.GetArcCostForVehicle(
            #     previous_index, index, vehicle_id)
            # tmp.append((manager.IndexToNode(previous_index),manager.IndexToNode(index)))
            tmp.append(previous_index)
            
            
        routes[vehicle_id] = np.array(tmp[1:])
        plan_output += ' {}\n'.format(manager.IndexToNode(index))
        # plan_output += 'Distance of the route: {}m\n'.format(route_distance)
        print(plan_output)
        # total_distance += route_distance
    # print('Total Distance of all routes: {}m'.format(total_distance))
    return routes

#
# parameters and reads
# while True:
len_avg_rds = 20000
sig_rds_len = 0

len_dna = 120000

overlap = 12000
sig_ovlp = 0

snp_rate = 0.01
num_ploid = 2
len_repeat = 0
num_repeats = 0

error_rate = 0.001

eos = True

max_error_n = round(error_rate*len_avg_rds*num_ploid*(1+(len_dna+num_repeats*len_repeat-len_avg_rds)/(len_avg_rds - overlap)))

data_path = '/mnt/d/projects/vrpassembler/sim_only/files'

paf_path = data_path + '/all_ovl.paf'

reads, rfa_file= generate_reads(len_avg_rds, sig_rds_len, len_dna, overlap, sig_ovlp, snp_rate, 
                                                    num_ploid, len_repeat, num_repeats, 
                                                    error_rate, data_path, eos)

avg_rl = 0
for i in reads:
    avg_rl += len(i)
avg_rl /= len(reads)
print(f'\naverage length of reads: {len_avg_rds}')
print(f'average length of overlap: {overlap}\n')

#
# run minimap2
process = subprocess.Popen(f"minimap2 -c -k 13 --for-only -D -t 6 \
                                -x ava-pb {rfa_file} {rfa_file} > {paf_path}",
                            shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
process.wait()
command_output = process.stdout.read().decode('utf-8')
print(command_output)


#%%
# initialize overlap mat
reload(overlap_helper_paf)

pen = 10
overlap_,tmp = overlap_helper_paf.process_paf(paf_path, len(reads), pen, avg_rl)
overlap = overlap_.copy()
#
ovmin = overlap.min()
m1 = np.identity(overlap.shape[0]-1)
if ovmin<0:
    overlap[1:,1:] -= ovmin
    
    m1 *= ovmin
    overlap[1:,1:] += m1
else:
    ovmax = overlap.max()
    tmp1 = overlap[1:,1:].copy()
    tmp2 = np.where(tmp1 == 0, ovmax, tmp1)
    m1 *= ovmax
    tmp2 -= m1
    overlap[1:,1:] = tmp2.copy()

dest_num = overlap.shape[0]


# run ortools
data = {}
data['distance_matrix'] = overlap

data['num_vehicles'] = num_ploid
data['depot'] = 0

manager = pywrapcp.RoutingIndexManager(overlap.shape[0], num_ploid, 0)

routing = pywrapcp.RoutingModel(manager)

def distance_callback(from_index, to_index):
    """Returns the distance between the two nodes."""
    # Convert from routing variable Index to distance matrix NodeIndex.
    from_node = manager.IndexToNode(from_index)
    to_node = manager.IndexToNode(to_index)
    return overlap[from_node,to_node]

transit_callback_index = routing.RegisterTransitCallback(distance_callback)

routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

search_parameters = pywrapcp.DefaultRoutingSearchParameters()
search_parameters.first_solution_strategy = (
    routing_enums_pb2.FirstSolutionStrategy.LOCAL_CHEAPEST_INSERTION)
search_parameters.local_search_metaheuristic = (
    routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH)
search_parameters.time_limit.seconds = 1
search_parameters.log_search = True

solution = routing.SolveWithParameters(search_parameters)

if solution:
    routes = print_solution(data, manager, routing, solution)
else:
    print('No solution found !')
    breakpoint()

#
# 判断结果，但不合并成序列，只看访问顺序

print('num ploids:',num_ploid)
print('num reads:',len(reads))
print('num vars:',num_ploid*len(reads)**num_ploid,'; min:',len(reads)**num_ploid)
print('error rate:',error_rate)
print('SNPs rate:',snp_rate)

fullmatch = True
for i in range(num_ploid):
    s1 = sum(routes[i][1:] - routes[i][:-1]-1)
    if s1 != 0:
        fullmatch = False
        break

print('SNPs recover:',fullmatch)

if not fullmatch:
    print('*'*50)
else:
    print('='*50)
    
    
