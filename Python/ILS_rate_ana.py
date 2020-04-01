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

best_ = filtered.loc[filtered['rate'] == 0.05555555555555555]['score']
second_best = filtered.loc[filtered['rate'] == 0.3888888888888888]['score']
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


    
