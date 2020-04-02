#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 10:58:05 2020

@author: timo
"""
import pandas as pd
import numpy as np
from scipy.stats import shapiro
from scipy.stats import mannwhitneyu, ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns


data = pd.read_csv('ILS_rate_exp.csv')

# =============================================================================
# FILTER AND SORT
# =============================================================================
filtered = data.drop(['Unnamed: 0', 'iter', 'time', 'passes'], axis=1)
outcome_measures = filtered.groupby(by=['iteration', 'rate']).agg(np.mean).reset_index()
filtered2 = outcome_measures.drop(['iteration'], axis=1)
outcome = filtered2.groupby(by='rate').agg([np.mean, np.std]).reset_index()



filtered3 = filtered2.loc[filtered2['rate'] != 0.0]
xticks = [round(x, 3) for x in filtered3['rate'].unique()]
# =============================================================================
# PLOT
# =============================================================================
plt.figure()
color = sns.color_palette('gray')
sns.boxplot(x='rate', y='score', data=filtered3, color='gray')
plt.xticks(np.arange(10), xticks, rotation=45)
plt.ylabel('Cutsize', fontsize=12)
plt.xlabel('Mutation rate', fontsize=12)
plt.ylim((58, 73))
plt.savefig('rateplot.png', dpi=300)
plt.show()

# =============================================================================
# GET BEST SCORES
# =============================================================================
scores = list(reversed(sorted(list(outcome[('score', 'mean')]))))
best_score = scores.pop()
second_best_score = scores.pop()

best_rate = float(outcome.loc[outcome[('score', 'mean')] == best_score]['rate'])
second_best_rate = float(outcome.loc[outcome[('score', 'mean')] == second_best_score]['rate'])

print(f'Best rates: {best_rate}, {second_best_rate}')

best_ = filtered3.loc[filtered3['rate'] == best_rate]['score']
second_best = filtered3.loc[filtered3['rate'] == second_best_rate]['score']
alphas = []
for i, data in enumerate([best_, second_best]):
    stat, p = shapiro(data)
    # print('Statistics=%.3f, p=%.3f' % (stat, p))
    # interpret
    alpha = 0.05
    if p > alpha:
     	print('Sample looks Gaussian (fail to reject H0)')
    else:
     	print('Sample does not look Gaussian (reject H0)')

print(ttest_ind(best_, second_best))
# statistic, p = mannwhitneyu(best_, second_best)


    
