#!/usr/bin/env pypy
# -*- coding: utf-8 -*-
from Graph import Graph
from Node import Node
from Bucket import Bucket
from Net import Net
import time
import random
from copy import deepcopy


def get_best_solution(graph):
    sol = [n.belongs_to_partition for n in graph.node_list]
    return sol

def create_solution(string_length):
    zeros = [0] * 250
    ones = [1] * 250
    z_o = zeros + ones

    sol = random.sample(z_o, string_length)
    
    return sol

        
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
    graph.compute_initial_gains()
    score = graph.calc_gain_sum()
#    print(score)
    
    best_gain_sum = score
    best_score_graph = graph
    
    while graph.partition_size > 0:
        node_to_change_A = graph.get_best_node(0)
        graph.flip_partition(node_to_change_A.index)
        
        node_to_change_B = graph.get_best_node(1)
        graph.flip_partition(node_to_change_B.index)
        
        graph.compute_gains(node_to_change_A)
        graph.compute_gains(node_to_change_B)
        
        score = graph.calc_gain_sum()
#        print(score)
        
        # Keep track of best gainSum and its accompanying graph
        if score > best_gain_sum:
            print('hoi')
            best_gain_sum = score
            best_score_graph = deepcopy(graph)
        
    best_score = best_score_graph.count_connections()
    best_score_solution = get_best_solution(best_score_graph)
    end_score = graph.count_connections()
    end_solution = get_best_solution(graph)
    
    return best_score, best_score_solution, end_score, end_solution

def local_search(graaf, stop):
#    best_graph = graaf
#    best_score = 2564
    scores = []
    best_scores = []
    
    for i in range(stop):
        score, solution, sc, sol = single_FM_run(graaf, 16)
#        print(sol[:5], solution[:5])
        graaf.define_partitions(sol)
        
        scores.append(sc)
        best_scores.append(score)
        
    return scores, best_scores
        

def MLS_test(node_list, iterations, stop):
    i = 0
    results = []

    graaf = Graph(node_list)
    solution = create_solution(500)
    graaf.define_partitions(solution)
    
    while i < iterations:
        result, result_ = local_search(graaf, stop)
        i += stop
        results.append(result)
        print(f'It {i}, score {result}, best_score{result_}')
        solution = create_solution(500)
        graaf.define_partitions(solution)
 

    
    
def MLS(node_list:list, maxcon:int, iterations:int):
    score = 2564

    i = 0
    result_dict = {key:[] for key in ['score', 'time']}
    
    start = time.time()
    graaf = Graph(node_list)
    solution = create_solution(500)
    
    localcounter = 0
    
    while i < iterations:
        start = time.time()
        
        graaf.define_partitions(solution)
        
        score_, solution_ = single_FM_run(graaf, maxcon)
        i += 1
        
        dur = time.time() - start
        result_dict['score'].append(score)
        result_dict['time'].append(dur)
        print(f"Iteration: {i}, Score: {score}, Time: {dur}s.")
        
        if score_ < score:
            localcounter += 1
            score, solution = score_, solution_  
#            print(f"Iteration: {i}, Score: {score}, Time: {dur}s, counter: {localcounter}")

        else:
            localcounter = 0
            # Create new solution, aka start anew
            start = time.time()
            graaf.define_partitions(create_solution(500))
            score, solution = single_FM_run(graaf, maxcon)
            i += 1
            
            dur = time.time() - start
            result_dict['score'].append(score)
            result_dict['time'].append(dur)
            print(f"Iteration: {i}, Score: {score}, Time: {dur}s.")
    
    return result_dict

def ILS(node_list:list, maxcon:int, iterations:int):
    score = 2564
    
    i=0
    
    start = time.time()
    graaf = Graph(node_list)
    solution = create_solution(500)
    
    

nodes, maxcon = parse_graph()

# Print nodes and connections
#for n in nodes:
#    print(f'{n.index}, {len(n.connection_locations)}, {n.connection_locations}')

#res = MLS(nodes, maxcon, 100)
MLS_test(nodes, 10, 10)

# Test of het iets te maken heeft met hergebruiken van solutions
#graaf = Graph(nodes)
#for i in range(300):
#    graaf.define_partitions(create_solution(500))
#    sc, so = single_FM_run(graaf, maxcon)
#    print(f'it {i}, score {sc}')
    

        
        
        
        
        
        
    
    
        
    
    