#!/usr/bin/env pypy
import search_algos
import numpy as np
from multiprocessing import Pool
from itertools import repeat
import pandas as pd
MAX_ITER = 10000
nodes = search_algos.parse_graph()

rates = list(np.linspace(0.0, 0.5, num=10))
rates = [float(x) for x in rates]

pool = Pool(10)

args = list(zip(repeat(nodes, 10), repeat(MAX_ITER, 10), rates))
ILS_iter_res12 = pool.starmap(search_algos.ILS_iter, args)
results = []
for i, dict_ in enumerate(ILS_iter_res12):
    dict_['iteration'] = repeat(i, len(dict_[list(dict_.keys())[0]]))
    dict_['rate'] = repeat(rates[i], len(dict_[list(dict_.keys())[0]]))
    results.append(pd.DataFrame(dict_))
results_total = pd.concat(results)
results_total.to_csv('ILS_rate_exp.csv')
    