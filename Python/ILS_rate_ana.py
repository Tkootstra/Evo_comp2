#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 10:58:05 2020

@author: timo
"""
import pandas as pd
import numpy as np
from scipy.stats import shapiro
from scipy.stats import mannwhitneyu


data = pd.read_csv('ILS_rate_exp.csv')

filtered = data.drop(['Unnamed: 0', 'iter', 'time', 'passes', 'iteration'], axis=1)
outcome_measures = filtered.groupby(by="rate").agg([np.mean, np.std]).reset_index()
scores = list(reversed(sorted(list(outcome_measures[('score', 'mean')]))))
best_score = scores.pop()
second_best_score = scores.pop()

best_rate = float(outcome_measures.loc[outcome_measures[('score', 'mean')] == best_score]['rate'])
second_best_rate = float(outcome_measures.loc[outcome_measures[('score', 'mean')] == second_best_score]['rate'])

best_ = filtered.loc[filtered['rate'] == best_rate]['score']
second_best = filtered.loc[filtered['rate'] == second_best_rate]['score']
alphas = []
for i, data in enumerate([best_, second_best]):
    stat, p = shapiro(data)
    print('Statistics=%.3f, p=%.3f' % (stat, p))
    # interpret
    alpha = 0.05
    if p > alpha:
     	print('Sample looks Gaussian (fail to reject H0)')
    else:
     	print('Sample does not look Gaussian (reject H0)')

statistic, p = mannwhitneyu(best_, second_best)


    
