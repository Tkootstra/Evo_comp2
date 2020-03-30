#!/usr/bin/env pypy
# -*- coding: utf-8 -*-
from Graph import Graph
from Node import Node
from Bucket import Bucket
from Net import Net
import time
import random
from copy import deepcopy
import numpy as np


def get_best_solution(node_list):
    sol = [n.belongs_to_partition for n in node_list]
    return sol

def create_solution(string_length):
    zeros = [0] * int((string_length /2))
    ones = [1] * int((string_length / 2))
    z_o = zeros + ones
    
    random.shuffle(z_o)
    
    return z_o

def mutate_solution(solution, rate):
    list0 = [i for i, v in enumerate(solution) if v == 0]
    list1 = [i for i, v in enumerate(solution) if v == 1]
    how_many = round(rate * len(solution))
    indices0 = list(np.random.choice(list0, size = how_many, replace=False))
    indices1 = list(np.random.choice(list1, size= how_many, replace=False))
    #endlist = [flip_bit(v) for i,v in enumerate(solution) if i in indices0 or i in indices1]
    
    endlist = []
    for i, v in enumerate(solution):
        if i in indices0 or i in indices1:
            endlist.append(flip_bit(v))
        else:
            endlist.append(v)

    return endlist
    
def flip_bit(bit):
    return 0 if bit==1 else 1

        
def parse_graph():
    maxcon = 0
    
    nodes = []
    file = 'Graph500.txt'
    txt = open(file, "r")
    lines = txt.readlines()
    for line in lines:
        conlocs = []
        splitted = line.split()
        index = int(splitted[0])-1
        ncon = int(splitted[2])
        [conlocs.append(int(x)-1) for x in splitted[3:len(splitted)]]
        new_node = Node(index, ncon, conlocs)
        nodes.append(new_node)
        
        # Counts max degree for bucket later on
        if ncon > maxcon:
            maxcon = ncon
    
    return nodes, maxcon
    

def single_FM_run(graph, maxcon):               
    # Compute initial gains 
    start = time.time()
    
    graph.compute_initial_gains()
    # print(f'gains {time.time() - start}')
    
    start = time.time()
    score = graph.calc_gain_sum()
    # print(f'gainsum {time.time() - start}')
    
#    best_gain_sum = score
#    best_score = None
#    best_score_solution = None
    

    start = time.time()
    while graph.free_nodes > 0:
        start0 = time.time()
        node_to_change_A = graph.get_best_node(0)
        getnode0 = time.time() - start0
        # print(f'getnode 0 {time.time() - start1}')
        
        start1 = time.time()
        graph.flip_partition(node_to_change_A.index)
        flip0 = time.time() - start1
        
        start1 = time.time()
        node_to_change_B = graph.get_best_node(1)
        getnode1 = time.time() - start1
        # print(f'getnode 1 {time.time() - start1}')
        
        start1 = time.time()
        graph.flip_partition(node_to_change_B.index)
        flip1 = time.time() - start1
        
        start1 = time.time()
        graph.compute_gains(node_to_change_A)
        gain0 = time.time() - start1
        
        start1 = time.time()
        graph.compute_gains(node_to_change_B)
        gain1 = time.time() - start1
        
                
        
        looptime = time.time() - start0
        if looptime > 1:
            print(f'{looptime}, gn {getnode0}, gn {getnode1}, f {flip0}, f {flip1}, g {gain0}, g {gain1}, \n')
        
#        print(f'loop {time.time() - start}')
       
#    start = time.time()
    end_score = graph.count_connections()
#    start = time.time()
    end_solution = get_best_solution(graph.node_list)
#     print(" getting best solution: {}".format(time.time() - start))
    
    
    return end_score, end_solution

def local_search(node_list, solution):
    begin = time.time()
    best = 9999 

    i = 0
    while(True):        
        graaf = Graph(node_list, solution)
        
        sc, solution = single_FM_run(graaf, 16)        
        i += 1
        
        if sc > best:
            break
        else:
            best = sc
#        print('single FM: {}'.format(time.time() - start))
    
    dur = time.time() - begin

    return best, i, dur, solution
        

def MLS(node_list, iterations):
    print('Running MLS...')
    
    result_dict = {key: [] for key in ['iter', 'score', 'time']}
    i = 0

    solution = create_solution(500)
    
    while i < iterations:
        result, i_, dur, sol = local_search(node_list, solution)
        i += i_

        solution = create_solution(500)

        result_dict['iter'].append(i)
        result_dict['score'].append(result)
        result_dict['time'].append(dur)
        print(f'{i}: {result}, {round(dur, 3)}')
    
    return result_dict

def ILS(node_list:list, iterations:int):
    print('Running ILS...')
    
    result_dict = {key: [] for key in ['iter', 'score', 'time']}
    start = time.time()

    i=0
    
    solution = create_solution(500)
    
    prev, i_, _, prev_solution = local_search(node_list, solution)
    temp = mutate_solution(prev_solution, 0.1)
    i+=i_
    dur = time.time() - start
    
    result_dict['iter'].append(i)
    result_dict['score'].append(prev)
    result_dict['time'].append(dur)
    
    while i < iterations:
        start = time.time()
        current, i_, _, current_solution = local_search(node_list, temp)
        i+=i_
        
        if current < prev:
            temp = current_solution
            prev = current
        else:
            temp = mutate_solution(current_solution, 0.1)

        dur = time.time() - start
        
        result_dict['iter'].append(i)
        result_dict['score'].append(current)
        result_dict['time'].append(dur)
        
        print(f'{i}: {current}, {round(dur, 3)}')
        # mutate je prev solution
        # Als dit beter is, ga verder met dit resultaat
        # Anders: mutate prev solution opnieuw en doe bovenstaand opnieuw
    
    return result_dict



nodes, maxcon = parse_graph()

#res_MLS = MLS(nodes, 100)
res_ILS = ILS(nodes, 1000)




    

        
        
        
        
        
        
    
    
        
    
    