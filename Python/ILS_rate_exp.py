#!/usr/bin/env pypy
import search_algos
import numpy as np
from multiprocessing import Pool
from itertools import repeat
import pandas as pd
MAX_ITER = 5000
nodes = search_algos.parse_graph()

rates = list(np.arange(0.0, 0.04, step=0.004))
rates = [float(x) for x in rates]
n_repeats = 5


pool = Pool(10)
total_res = []
args = list(zip(repeat(nodes, 10), repeat(MAX_ITER, 10), rates))

for x in range(n_repeats):
    ILS_iter_res12 = pool.starmap(search_algos.ILS_iter, args)
    results = []
    for i, dict_ in enumerate(ILS_iter_res12):
        dict_['iteration'] = repeat(x, len(dict_[list(dict_.keys())[0]]))
        dict_['rate'] = repeat(rates[i], len(dict_[list(dict_.keys())[0]]))
        results.append(pd.DataFrame(dict_))
    results_total = pd.concat(results)
    total_res.append(results_total)
all_ = pd.concat(total_res)
all_.to_csv('ILS_rate_exp.csv')
    
    

    