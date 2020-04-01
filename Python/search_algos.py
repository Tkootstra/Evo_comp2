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


def get_solution_from_nodes(node_list):
    sol = [n.belongs_to_partition for n in node_list]
    return sol

def create_solution(string_length):
    zeros = [0] * int((string_length /2))
    ones = [1] * int((string_length / 2))
    z_o = zeros + ones
    
    random.shuffle(z_o)
    
    return z_o

def hamming_distance(sol1, sol2):
    distance = 0
    for bit1, bit2 in zip(sol1, sol2):
        if bit1 != bit2:
            distance+=1
    return distance

def invert_solution(solution):
    return list(map(flip_bit, solution))

def uniform_crossover(parent1, parent2):
    
    if hamming_distance(parent1, parent2) > (len(parent1)/2):
        parent1 = invert_solution(parent1)
    child = [None] * len(parent1)
    sample_indices = []
    samplepool = []
    for i, (bit_a, bit_b) in enumerate(zip(parent1, parent2)):
        if bit_a == bit_b:
            child[i] = bit_a
        else:
            sample_indices.append(i)
    # compute balance
    zero_amount = len(list(filter(lambda x: True if x == 0 else False, child)))
    one_amount = len(list(filter(lambda x: True if x == 1 else False, child)))
    amount_of_zeros = int((len(parent1) / 2) - one_amount)
    amount_of_ones = int((len(parent1) / 2) - zero_amount)
    ones = [1] * amount_of_ones
    zeroes = [0] * amount_of_zeros
    
    sample_pool = ones + zeroes
    random.shuffle(sample_pool)
    
    for missing in sample_indices:
        child[missing] = sample_pool.pop()
    return child

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
        
    
    return nodes
    

def single_FM_run(graph): 
    best_score = 9999
    best_gainsum = -9999
    best_sol = None
              
    # Compute initial gains     
    graph.compute_initial_gains()

    start = time.time()
    while graph.free_nodes > 0:
        start0 = time.time()
        node_to_change_A = graph.get_best_node(0)
#        getnode0 = time.time() - start0
        
#        start1 = time.time()
        graph.flip_partition(node_to_change_A.index)
#        flip0 = time.time() - start1
        
#        start1 = time.time()
        node_to_change_B = graph.get_best_node(1)
#        getnode1 = time.time() - start1
        
#        start1 = time.time()
        graph.flip_partition(node_to_change_B.index)
#        flip1 = time.time() - start1
        
#        start1 = time.time()
        graph.compute_gains(node_to_change_A)
#        gain0 = time.time() - start1
        
#        start1 = time.time()
        graph.compute_gains(node_to_change_B)
#        gain1 = time.time() - start1
        
#        start1 = time.time()
        score = graph.count_connections()
                
        if score < best_score:            
            best_score = score
            best_sol = get_solution_from_nodes(graph.node_list)
        
#        counttime = time.time() - start1
        
#        looptime = time.time() - start0
#        if graph.free_nodes % 100 == 0:
#            print(f'loop {looptime}') #: ct {counttime}, gn {getnode0}, gn {getnode1}, f {flip0}, f {flip1}, g {gain0}, g {gain1}, \n')
              
    return best_score, best_sol

def local_search(node_list, solution, i, iterations):
    stop = False
    if i >= (iterations):
        return None, i, None, None
    begin = time.time()
    best = 9999 

    local_i = 0
    
    while(i < iterations):     
        start = time.time()
        graaf = Graph(node_list, solution)
        
        sc, solution = single_FM_run(graaf)        
        local_i += 1
        i += 1
        
        if sc >= best: # score is worse
            break
        else:
            best = sc
#        print('single FM: {}, score {}'.format(time.time() - start, sc))
    
    dur = time.time() - begin

    return best, i, dur, solution 

      

def MLS_iter(node_list, iterations):
    print('Running MLS...')
    
    result_dict = {key: [] for key in ['iter', 'score', 'time', 'passes']}

    i = 0
    
    solution = create_solution(500)
    
    while i < iterations:
        prev_i = i
        result, i_, dur, sol = local_search(node_list, solution, i, iterations)
        i = i_
        delta_i = i - prev_i
        solution = create_solution(500)

        result_dict['iter'].append(i)
        result_dict['score'].append(result)
        result_dict['time'].append(dur)
        result_dict['passes'].append(delta_i)
        print(f'{i}: Endscore: {result}, Mean duration: {round(dur / i_, 3)}')
    
    return result_dict


def MLS_time(node_list, iterations, stop_time):
    print('Running MLS...')
    
    result_dict = {key: [] for key in ['iter', 'score', 'time', 'passes']}

    i = 0
    total_dur = 0
    print(total_dur)
    solution = create_solution(500)
    
    while i < iterations and total_dur < stop_time:
        prev_i = i
        result, i_, dur, sol = local_search(node_list, solution, i, iterations)
        i = i_
        delta_i = i - prev_i
        
        solution = create_solution(500)

        result_dict['iter'].append(i)
        result_dict['score'].append(result)
        result_dict['time'].append(dur)
        result_dict['passes'].append(delta_i)
        # print(f'{i}: Endscore: {result}, Mean duration: {round(dur / i_, 3)}')
        total_dur += dur

    
    return result_dict



def ILS_iter(node_list:list, iterations:int, rate=0.1):
    print('Running ILS...')
    
    result_dict = {key: [] for key in ['iter', 'score', 'time', 'passes']}
    start = time.time()

    i = 0
    
    solution = create_solution(500)
    prev_i = i
    prev, i_, _, prev_solution = local_search(node_list, solution, i, iterations)
    temp = mutate_solution(prev_solution, rate)
    i = i_
    delta_i = i - prev_i
    dur = time.time() - start
    
    result_dict['iter'].append(i)
    result_dict['score'].append(prev)
    result_dict['time'].append(dur)
    result_dict['passes'].append(delta_i)
    
    # print(f'{i}: Endscore: {prev}, Mean duration: {round(dur / delta_i, 3)} ({delta_i} passes)')
    
    while i < iterations:
        start = time.time()
        prev_i = i
        current, i_, _, current_solution = local_search(node_list, temp, i, iterations)
        i = i_
        delta_i = i - prev_i
        
        if current < prev:
            temp = current_solution
            prev = current
        else:
            temp = mutate_solution(temp, rate)

        dur = time.time() - start
        
        result_dict['iter'].append(i)
        result_dict['score'].append(current)
        result_dict['time'].append(dur)
        result_dict['passes'].append(delta_i)
        
        
        # print(f'{i}: Endscore: {current}, Mean duration: {round(dur / delta_i, 3)} ({delta_i} passes)')
        # mutate je prev solution
        # Als dit beter is, ga verder met dit resultaat
        # Anders: mutate prev solution opnieuw en doe bovenstaand opnieuw
    
    return result_dict

def ILS_time(node_list:list, iterations:int, stop_time:int, rate=0.1):
    print('Running ILS...')
    
    result_dict = {key: [] for key in ['iter', 'score', 'time', 'passes']}
    start = time.time()

    i = 0
    total_dur = 0
    solution = create_solution(500)
    prev_i = i
    prev, i_, _, prev_solution = local_search(node_list, solution, i, iterations)
    temp = mutate_solution(prev_solution, rate)
    i = i_
    delta_i = i - prev_i
    dur = time.time() - start
    
    
    
    print(f'{i}: Endscore: {prev}, Mean duration: {round(dur / delta_i, 3)} ({delta_i} passes)')
    result_dict['iter'].append(i)
    result_dict['score'].append(prev)
    result_dict['time'].append(dur)
    result_dict['passes'].append(delta_i)
    
    while i < iterations and total_dur < stop_time:
        start = time.time()
        prev_i = i
        current, i_, _, current_solution = local_search(node_list, temp, i, iterations)
        
        i = i_
        delta_i = i - prev_i
        
        if current < prev:
            temp = current_solution
            prev = current
        else:
            temp = mutate_solution(temp, rate)

        dur = time.time() - start
        
        result_dict['iter'].append(i)
        result_dict['score'].append(current)
        result_dict['time'].append(dur)
        result_dict['passes'].append(delta_i)
        global_dur = time.time() - start
       

        print(f'{i}: Endscore: {current}, Mean duration: {round(dur / delta_i, 3)} ({delta_i} passes)')
        total_dur += global_dur
        # mutate je prev solution
        # Als dit beter is, ga verder met dit resultaat
        # Anders: mutate prev solution opnieuw en doe bovenstaand opnieuw
    
    return result_dict


def GLS_iter(node_list:list, iterations:int, population_size=50):
    print('Running GLS...')
    result_dict = {key: [] for key in ['iter', 'score', 'time', 'event','passes']}
    population = []
    for ii in range(population_size):
        population.append(create_solution(500))
    
    i = 0
    cont = True
    while cont:
        begin = time.time()
        worst_score = 0
        
        for idx, solution in enumerate(population): 
            # check the performance of all solution first
            prev_i = i
            result, i_, dur, sol = local_search(node_list, solution, i, iterations)
            if result is None:
                cont = False
                print('breaking')
                break
            
            
            i = i_  
            delta_i = i - prev_i
            result_dict['iter'].append(i)
            result_dict['score'].append(result)
            result_dict['time'].append(dur)
            result_dict['event'].append('normal')
            result_dict['passes'].append(delta_i)
            if result > worst_score:
                worst_score = result
                worst_idx = idx
            
        # check worst score and uniform x over child
        # if child is better, replace worst with child
            
        if i < iterations:
            begin = time.time()
            pa1 = random.randint(0, population_size-1)
            pa2 = random.randint(0, population_size-1)
            while pa1 == pa2:
                pa2 = random.randint(0, population_size-1)
            child = uniform_crossover(population[pa1], population[pa2])
            prev_i = i
            child_score, i_, dur, _  = local_search(node_list, child, i, iterations)
            i = i_
            delta_i = i - prev_i
            if child_score <= worst_score:
                population[worst_idx] = child
            dur = time.time() - begin
            result_dict['iter'].append(i)
            result_dict['score'].append(child_score)
            result_dict['time'].append(dur)
            result_dict['event'].append('crossover')
            result_dict['passes'].append(delta_i)
    return result_dict
    
def GLS_time(node_list:list, iterations:int, stop_time, population_size=50):
    print('Running GLS...')
    result_dict = {key: [] for key in ['iter', 'score', 'time', 'event','passes']}
    population = []
    for ii in range(population_size):
        population.append(create_solution(500))
    total_dur = 0
    i = 0
    cont = True
    while cont:
        begin = time.time()
        worst_score = 0
        
        for idx, solution in enumerate(population): 
            # check the performance of all solution first
            prev_i = i
            result, i_, dur, sol = local_search(node_list, solution, i, iterations)
            
            
            
            i = i_  
            delta_i = i - prev_i
            result_dict['iter'].append(i)
            result_dict['score'].append(result)
            result_dict['time'].append(dur)
            result_dict['event'].append('normal')
            result_dict['passes'].append(delta_i)

            if result > worst_score:
                worst_score = result
                worst_idx = idx
            total_dur += dur
            if total_dur >= stop_time:
                cont = False
                print('breaking')
                break

        # check worst score and uniform x over child
        # if child is better, replace worst with child
            
        if i < iterations and total_dur < stop_time:
            begin = time.time()
            pa1 = random.randint(0, population_size-1)
            pa2 = random.randint(0, population_size-1)
            while pa1 == pa2:
                pa2 = random.randint(0, population_size-1)
            child = uniform_crossover(population[pa1], population[pa2])
            prev_i = i
            child_score, i_, dur, _  = local_search(node_list, child, i, iterations)
            i = i_
            delta_i = i - prev_i
            if child_score <= worst_score:
                population[worst_idx] = child
            dur = time.time() - begin
            result_dict['iter'].append(i)
            result_dict['score'].append(child_score)
            result_dict['time'].append(dur)
            result_dict['event'].append('crossover')
            result_dict['passes'].append(delta_i)
            total_dur += time.time() - begin
    return result_dict



# res_MLS = MLS_time(nodes, 99999999, 20)
# res_ILS = ILS_time(nodes, 999999999, 20, rate=0.1)
# gls = GLS_time(nodes, 99999999, 50, 20)




 
    
    
        
    
    