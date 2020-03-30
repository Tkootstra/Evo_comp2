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

    # sol = random.sample(z_o, string_length)
    
    random.shuffle(z_o)
    
    return z_o

        
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
    


# def single_swap_FM(eerste_random_solution, graph):
#     # 1. maak aparte lijsten voor partitie a en b
#     # 2. voor elke node: check of die fixed is, voeg dan pas toe
#     # 3. calculate base gain voor de hele graph/node list
#     # 3. voor elke partitie lijst, doe:
         
#         - check elke node: doe de swap, calculate score <- niet nodig?
#         - nieuwe gain = calculate gain
#         - best node = kies een random beste gain
#         - maak de swap voor de beste gain node
#         - fix deze node
        
#         - return de graph die eruit komt
        
# def FM_pass(solution, swapcounter):
#     - maak alle nodes free
    
#     - zolang er vrije nodes zijn en swapcounter < 10,000:
#         - graph = single_swap_FM(solution, graph)
#         - swapcounter += 1
        
#         - bereken score
#         - houd beste score bij
        
#     - return best_score, swapcounter
    
# def FM(graph, swapcounter):
#     - score = bereken score
#     - graph, swapcounter = FM_pass(graph, swapcounter)
#     - score2 = bereken score
    
#     - while(score2 < score):  # zolang het beter wordt
#         - score = score2
#         - graph, swapcounter = FM_pass(graph, swapcounter)
#         - score2 = bereken score
        



def single_FM_run(graph, maxcon):       
    
    # # zoals yannick
    # 1. FM swap (args = solution, fixed_list)
    #     - maak aparte lijsten van de parttities
    #     - voeg alleen een node toe aan de lijst als hij niet fixed is
    #     Voor Zowel partitie A als partitie B, doe:
    #         - calculate pre score (base gain)
    #         - doe meteen een swap op deze bit
    #         - je krijgt een nieuwe solution terug
    #         - calculate de post score (gain)
    #         - prescore - postscore = gain
    #         - hou de gain bij voor elke swap
    #         - kies een random beste gain (als er meerdere zijn, anders return je gewoon de beste)
    #         - maak de swap op de solution die je gekozen hebt
    #         - de node die je geswapt heb wordt nu fixed
    #     - return solution die terugkomt uit de swap (alleen voor partitie B?????)
        
    # 2. FM pass(solution, swapcounter):
    #     - maak eerst alle nodes free
    #     - zolang er nog vrije nodes zijn en swapcounter < 10000:
    #         solution, fixed_list = fm_swap(solution, fixed_list)
    #         swapcounter++
    #         score (gain) = getscore(solution)
    #         hou beste score bij
            
    #     return best_score, swapcounter
    
    # 3. FM(solution, swapcounter):
    #     score = getscore(solution)
    #     solution, swapcounter = FM pass(solution, swapcounter)
    #     while (getscore(solution) < score):
    #         score = getscore(solution)
    #         individual, swapcounter = FM pass(solution, swapcounter)
    #     else:
    #         log de score, en de swapcounter
            
    #     return solution, swapcounter
    
        
    # Compute initial gains 
    start = time.time()
    
    graph.compute_initial_gains()
    # print(f'gains {time.time() - start}')
    
    start = time.time()
    score = graph.calc_gain_sum()
    # print(f'gainsum {time.time() - start}')
#    print(score)
    
    best_gain_sum = score
    best_score = None
    best_score_solution = None
    
    
    
    
    
    start = time.time()
    while graph.free_nodes > 0:
        # print(f'free nodes {graph.free_nodes}')
        start1 = time.time()
        node_to_change_A = graph.get_best_node(0)
        # print("best node is {}".format(node_to_change_A.index))
        # print(f'getnode 0 {time.time() - start1}')
        graph.flip_partition(node_to_change_A.index)
        
        start1 = time.time()
        node_to_change_B = graph.get_best_node(1)
        # print(f'getnode 1 {time.time() - start1}')
        graph.flip_partition(node_to_change_B.index)
        # start = time.time()
        graph.compute_gains(node_to_change_A)
        graph.compute_gains(node_to_change_B)
        
        score = graph.calc_gain_sum()
        # print('compute gains: {}'.format(time.time() - start))
#        print(score)
        
        # Keep track of best gainSum and its accompanying graph
        
        
    # print(f'loops {time.time() - start1}')
       
#    best_score = best_score_graph.count_connections()
#    best_score_solution = get_best_solution(best_score_graph)
    start = time.time()
    end_score = graph.count_connections()
    start = time.time()
    end_solution = get_best_solution(graph.node_list)
    # print(" getting best solution: {}".format(time.time() - start))
    
    
    
    
    
    return end_score, end_solution

def local_search(node_list, solution):
    begin = time.time()
    best = 9999 
    # sc = 9998
    scores = []

    # print(solution[:10])
    i = 0
    while(True):
        start = time.time()
        
        graaf = Graph(node_list, solution)
        # print(f'graaf {time.time() - start}')
        
        sc, solution = single_FM_run(graaf, 16)
        # print('fblsfpibeflbflffbel')
        # print(f'{i}: {sc}, {time.time() - start}, {score}, {solution[:10]}')
        
        scores.append(sc)
        # best_scores.append(score)
        i += 1
        
        if sc > best:
            break
        else:
            best = sc
        print('single FM: {}'.format(time.time() - start))
    dur = time.time() - begin
    # print('done with current pass')    
    return best, i, dur, solution
        

def MLS(node_list, iterations):
    result_dict = {key: [] for key in ['iter', 'score', 'time']}
    i = 0
    results = []

    solution = create_solution(500)
    
    while i < iterations:
        result, i_, dur, sol = local_search(node_list, solution)
        i += i_
        results.append(result)
        # print(f'It {i}, score {result}, best_score{result_}')
        solution = create_solution(500)
#        graaf.define_partitions(solution)
        result_dict['iter'].append(i)
        result_dict['score'].append(result)
        result_dict['time'].append(dur)
        print(i, result, dur)
    return result_dict

def ILS(node_list:list, iterations:int):
    result_dict = {key: [] for key in ['iter', 'score', 'time']}
    begin = time.time()
    score = 2564
    best_score = None
    i=0
    begin = time.time()
    
    start = time.time()
    solution = create_solution(500)
    
    prev, i_, dur, prev_solution = local_search(node_list, solution)
    temp = mutate_solution(prev_solution, 0.1)
    best_score = prev
    i+=i_
    end = time.time()
    result_dict['iter'].append(i)
    result_dict['score'].append(prev)
    result_dict['time'].append(end-begin)
    
    while i < iterations:
        begin = time.time()
        current, i_, dur, current_solution = local_search(node_list, temp)
        if current < prev:
            temp = current_solution
            prev = current
            kek = "doing normal"
        else:
            temp = mutate_solution(current_solution, 0.1)
            kek = "doing mutation"
        i+=i_
        end = time.time()
        result_dict['iter'].append(i)
        result_dict['score'].append(current)
        result_dict['time'].append(end-begin)
        
        print(i, current, end-begin, kek)
        # mutate je prev solution
        # Als dit beter is, ga verder met dit resultaat
        # Anders: mutate prev solution opnieuw en doe bovenstaand opnieuw
    return result_dict
        
        
    
    
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



nodes, maxcon = parse_graph()

# res = ILS(nodes, 100)
res = MLS(nodes, 100)
# Print nodes and connections
#for n in nodes:
#    print(f'{n.index}, {len(n.connection_locations)}, {n.connection_locations}')

#res = MLS(nodes, maxcon, 100)
# res = MLS(nodes, 100)

# Test of het iets te maken heeft met hergebruiken van solutions
#graaf = Graph(nodes)
#for i in range(300):
#    graaf.define_partitions(create_solution(500))
#    sc, so = single_FM_run(graaf, maxcon)
#    print(f'it {i}, score {sc}')
    

        
        
        
        
        
        
    
    
        
    
    