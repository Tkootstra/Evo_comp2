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
MAX_TIME = 60
RATE = 0.1  # to be determined
MAX_ITER = 999999999

# # MLS: time
# pool = Pool(12)
# max_iters = [MAX_TIME] * 12
# # copies =  repeat(nodes.copy())
# args = list(zip(repeat(nodes, 12), repeat(MAX_ITER, 12), repeat(MAX_TIME, 12)))
# MLS_iter_res12 = pool.starmap(search_algos.MLS_time, args)
# pool = Pool(13)
# args = list(zip(repeat(nodes, 13), repeat(MAX_ITER, 13), repeat(MAX_TIME, 12)))
# MLS_iter_res13 = pool.starmap(search_algos.MLS_time, args)
# total = MLS_iter_res12 + MLS_iter_res13
# results = []
# for i, dict_ in enumerate(total):
#     dict_['iteration'] = repeat(i, len(dict_[list(dict_.keys())[0]]))
#     results.append(pd.DataFrame(dict_))
# results_total = pd.concat(results)
# results_total.to_csv('MLS_time.csv')

# del results_total, MLS_iter_res12, MLS_iter_res13, total, results

# # ILS time
# pool = Pool(12)
# max_iters = [MAX_ITER] * 12
# # copies =  repeat(nodes.copy())
# args = list(zip(repeat(nodes, 12), repeat(MAX_ITER, 12), repeat(MAX_TIME, 12)))
# ILS_iter_res12 = pool.starmap(search_algos.ILS_time, args)
# args = list(zip(repeat(nodes, 13), repeat(MAX_ITER, 13), repeat(MAX_TIME, 12)))
# ILS_iter_res13 = pool.starmap(search_algos.ILS_time, args)
# total = ILS_iter_res12 + ILS_iter_res13
# results = []
# for i, dict_ in enumerate(total):
#     dict_['iteration'] = repeat(i, len(dict_[list(dict_.keys())[0]]))
#     results.append(pd.DataFrame(dict_))
# results_total = pd.concat(results)
# results_total.to_csv('ILS_time.csv')

# del results_total, ILS_iter_res12, ILS_iter_res13, total, results


# GLS time
pool = Pool(12)
max_iters = [MAX_ITER] * 12
# copies =  repeat(nodes.copy())
args = list(zip(repeat(nodes, 12), repeat(MAX_ITER, 12), repeat(MAX_TIME, 12)))
GLS_iter_res12 = pool.starmap(search_algos.GLS_time, args)
args = list(zip(repeat(nodes, 13), repeat(MAX_ITER, 13)))
GLS_iter_res13 = pool.starmap(search_algos.GLS_time, args)
total = GLS_iter_res12 + GLS_iter_res13
results = []
for i, dict_ in enumerate(total):
    dict_['iteration'] = repeat(i, len(dict_[list(dict_.keys())[0]]))
    results.append(pd.DataFrame(dict_))
results_total = pd.concat(results)
results_total.to_csv('GLS_time.csv')

del results_total, GLS_iter_res12, GLS_iter_res13, total, results













