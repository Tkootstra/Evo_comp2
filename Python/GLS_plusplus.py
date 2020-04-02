#!/usr/bin/env pypy
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 15:34:25 2020

@author: timo
"""
import search_algos
from multiprocessing import Pool
from itertools import repeat
import pandas as pd

nodes = search_algos.parse_graph()
# PARAMS
# iterative experiments
MAX_ITER = 10000
RATE = 0.004  # to be determined
N_CORES = 12

# GLS++ iterative
pool = Pool(12)
max_iters = [MAX_ITER] * N_CORES
# copies =  repeat(nodes.copy())
args = list(zip(repeat(nodes, 12), repeat(MAX_ITER, 12)))
MLS_iter_res12 = pool.starmap(search_algos.GLS_own_iter, args)
pool = Pool(13)
args = list(zip(repeat(nodes, 13), repeat(MAX_ITER, 13)))
MLS_iter_res13 = pool.starmap(search_algos.GLS_own_iter, args)
total = MLS_iter_res12 + MLS_iter_res13
results = []
i = 0
for dict_ in total:
    dict_['iteration'] = repeat(i, len(dict_[list(dict_.keys())[0]]))
    results.append(pd.DataFrame(dict_))
    i+=1
results_total = pd.concat(results)
results_total.to_csv('GLS_plusplus.csv')

del results_total, MLS_iter_res12, MLS_iter_res13, total, results
