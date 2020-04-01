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
MAX_ITER = 50
RATE = 0.0555556  # to be determined
N_CORES = 12
# MLS: iterative
pool = Pool(12)
max_iters = [MAX_ITER] * N_CORES
# copies =  repeat(nodes.copy())
args = list(zip(repeat(nodes, 12), repeat(MAX_ITER, 12)))
MLS_iter_res12 = pool.starmap(search_algos.MLS_iter, args)
pool = Pool(13)
args = list(zip(repeat(nodes, 13), repeat(MAX_ITER, 13)))
MLS_iter_res13 = pool.starmap(search_algos.MLS_iter, args)
total = MLS_iter_res12 + MLS_iter_res13
results = []
for i, dict_ in enumerate(total):
    dict_['iteration'] = repeat(i, len(dict_[list(dict_.keys())[0]]))
    results.append(pd.DataFrame(dict_))
results_total = pd.concat(results)
results_total.to_csv('MLS_iter.csv')

del results_total, MLS_iter_res12, MLS_iter_res13, total, results

# ILS iter
pool = Pool(12)
max_iters = [MAX_ITER] * 12
# copies =  repeat(nodes.copy())
args = list(zip(repeat(nodes, 12), repeat(MAX_ITER, 12), repeat(RATE, 12)))
ILS_iter_res12 = pool.starmap(search_algos.ILS_iter, args)
pool = Pool(13)
args = list(zip(repeat(nodes, 13), repeat(MAX_ITER, 13), repeat(RATE, 12)))
ILS_iter_res13 = pool.starmap(search_algos.ILS_iter, args)
total = ILS_iter_res12 + ILS_iter_res13
results = []
for i, dict_ in enumerate(total):
    dict_['iteration'] = repeat(i, len(dict_[list(dict_.keys())[0]]))
    results.append(pd.DataFrame(dict_))
results_total = pd.concat(results)
results_total.to_csv('ILS_iter.csv')

del results_total, ILS_iter_res12, ILS_iter_res13, total, results


# GLS iter
pool = Pool(12)
max_iters = [MAX_ITER] * 12
# copies =  repeat(nodes.copy())
args = list(zip(repeat(nodes, 12), repeat(MAX_ITER, 12)))
GLS_iter_res12 = pool.starmap(search_algos.GLS_iter, args)
pool = Pool(13)
args = list(zip(repeat(nodes, 13), repeat(MAX_ITER, 13)))
GLS_iter_res13 = pool.starmap(search_algos.GLS_iter, args)
total = GLS_iter_res12 + GLS_iter_res13
results = []
for i, dict_ in enumerate(total):
    dict_['iteration'] = repeat(i, len(dict_[list(dict_.keys())[0]]))
    results.append(pd.DataFrame(dict_))
results_total = pd.concat(results)
results_total.to_csv('GLS_iter.csv')

del results_total, GLS_iter_res12, GLS_iter_res13, total, results













