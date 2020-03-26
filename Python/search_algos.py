#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from Graph import Graph
from Node import Node
from Bucket import Bucket
import time
import random


def init_gain(graph, current_bucket):
    for i in range(len(graph.node_list)):
        if not graph.node_list[i].is_fixed:
            current = graph.node_list[i]
            graph.node_list[current.index].flip_partition()
            
            gain = graph.count_single_cell_connections(current.index, 0)
            # print(gain)
            graph.node_list[current.index].flip_partition()
            
            current_bucket.add_to_bucket(current.belongs_to_partition, gain, current)
            graph.node_list[i].gain = gain
    return current_bucket
                 

def update_gain(graph, current_bucket, neighbors):
    for index in neighbors:
        if not graph.node_list[index].is_fixed:
            current = graph.node_list[index]
            graph.node_list[current.index].flip_partition()
            gain = graph.count_single_cell_connections(current.index, current.gain)
            graph.node_list[current.index].flip_partition()
            
            current_bucket.update_bucket(current.belongs_to_partition, gain, current)
            graph.node_list[index].gain = gain
    
def get_best_soltution(graph):
    sol = []
    for node in graph.node_list:
        sol.append(node.belongs_to_partition)
    return sol
            

def parse_graph():
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
    
    
    return nodes

def single_FM_run(graph):
    bucketA = dict()
    bucketB = dict()
    
    results = Bucket(bucketA, bucketB)
    
    bucket_A_size = results.bucketA_size
    
    results = init_gain(graph, results)
    score = results.gain_sum()
    
    best_gain_sum = score
    best_score_graph = graph
    
    while bucket_A_size > 0:
        node_to_change_A = results.pop_from_bucket_key(0)
        graph.node_list[node_to_change_A.index].flip_partition()
        
        node_to_change_B = results.pop_from_bucket_key(1)
        graph.node_list[node_to_change_B.index].flip_partition()
        results = update_gain(graph, results, node_to_change_A.connection_locations)
        results = update_gain(graph, results, node_to_change_B.connection_locations)
        score = results.gain_sum()
        bucket_A_size = results.bucketA
        
        if score > best_gain_sum:
            best_gain_sum = score
            best_score_graph = graph
        
    best_score = best_score_graph.count_connections()
    best_score_solution = get_best_soltution(best_score_graph)
    return best_score, best_score_solution

def create_solution(string_length):
    
    values = []
    for s in range(string_length):
        values.append(random.randint(0, 1))
#        values = np.random.randint(2, size=string_length)
    return values

    
def MLS(node_list:list, iterations:int):
    start = time.time()
    result_dict = {key:[] for key in ['best_score', 'time']}
    
    for i in range(iterations):
        graaf = Graph(node_list)
        graaf.define_partitions(create_solution(500))
        best_score, best_score_solution = single_FM_run(graaf)
        dur = time.time() - start
        result_dict['best_score'].append(best_score)
        result_dict['time'].append(dur)
        print("Iteration: {}, Score: {}, Time: {}, s.".format(i+1, best_score, dur))
    
    return result_dict

nodes = parse_graph()

res = MLS(nodes, 10)

        
        
        
        
        
        
    
    
        
    
    