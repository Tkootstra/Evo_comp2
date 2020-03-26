#!/usr/bin/env pypy
# -*- coding: utf-8 -*-
from Graph import Graph
from Node import Node
from Bucket import Bucket
from Net import Net
import time
import random


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

#def init_gain(graph, current_bucket):
#    # This works fine
#    for i in range(len(graph.node_list)):
#        current = graph.node_list[i]
#        
#        if not current.is_fixed:
#            graph.node_list[current.index].flip_partition()
#            gain = graph.count_single_cell_connections(current.index, current.gain)
#            graph.node_list[current.index].flip_partition()
#            
#            if gain < -16 or gain > 16:
#                print(gain)
#            
#            current_bucket.add_to_bucket(current.belongs_to_partition, gain, current)
#            graph.node_list[i].gain = gain
#    
#    return current_bucket
                  

#def update_gain(graph, current_bucket, neighbors):
#    # Gains get out of hand, usually towards -20 or -30. -16 should be minimum
#    for idx in neighbors:
#        current = graph.node_list[idx]
#        
#        if not current.is_fixed:
#            graph.node_list[current.index].flip_partition()                         # Temp flip
#            gain = graph.count_single_cell_connections(current.index, current.gain)
#            graph.node_list[current.index].flip_partition()                         # Flip back
#            
#            if gain < -16 or gain > 16:
#                print(f'{idx}: {current.gain} -> {gain} = {gain - current.gain}')
#            
#            current_bucket.update_bucket(current.belongs_to_partition, gain, current)
#            graph.node_list[idx].gain = gain
#            
#    return current_bucket
    


def single_FM_run(graph, maxcon):        
    # Create a new bucket
    results = Bucket(maxcon)
    
    # Compute initial gains    
#    results.init_gain(graph)
    results.init_gain_test(graph)
    score = results.gain_sum()
    
    best_gain_sum = score
    best_score_graph = graph
    
    while results.bucketA_size > 0:
        node_to_change_A = results.pop_from_bucket_key(0)
        graph.node_list[node_to_change_A.index].flip_partition()
        
        node_to_change_B = results.pop_from_bucket_key(1)
        graph.node_list[node_to_change_B.index].flip_partition()
        
#        results.update_gain(graph, node_to_change_A.connection_locations)
#        results.update_gain(graph, node_to_change_B.connection_locations)
        results.update_gain_test(graph, node_to_change_A.index)
        results.update_gain_test(graph, node_to_change_B.index)
        
        score = results.gain_sum()
        
        # Keep track of best gainSum and its accompanying graph
        if score > best_gain_sum:
            best_gain_sum = score
            best_score_graph = graph
        
    best_score = best_score_graph.count_connections()
    best_score_solution = get_best_solution(best_score_graph)
    
    return best_score, best_score_solution

    
def MLS(node_list:list, maxcon:int, iterations:int):
    score = 2564
    
    
    i = 0
    result_dict = {key:[] for key in ['score', 'time']}
    
    start = time.time()
    graaf = Graph(node_list)
    solution = create_solution(500)
    
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
            score, solution = score_, solution_  

        else:
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

nodes, maxcon = parse_graph()

# Print nodes and connections
#for n in nodes:
#    print(f'{n.index}, {len(n.connection_locations)}, {n.connection_locations}')

#res = MLS(nodes, maxcon, 10)

# Test of het iets te maken heeft met hergebruiken van solutions
graaf = Graph(nodes)
for i in range(300):
    graaf.define_partitions(create_solution(500))
    sc, so = single_FM_run(graaf, maxcon)
    print(f'it {i}, score {sc}')
    

        
        
        
        
        
        
    
    
        
    
    