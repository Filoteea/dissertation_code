# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 03:14:26 2022

@author: filot
"""

# plots models for 2018

# Standard Library imports
import argparse
import gzip
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import os
import pandas as pd
import sys
import xarray as xr
import csv
import random
import matplotlib.cm as cm
import scipy.linalg
import scipy.stats
from scipy.stats import spearmanr
from scipy.stats import pearsonr



# Third party imports
from collections import OrderedDict
from datetime import datetime
from numpy import genfromtxt


# Semi-local imports
import name_qch4_couple.io
import name_qch4_couple.name
import name_qch4_couple.plot_h2


# Local imports
import routines
import chem_co
import chem_ch4_validation
import merge_csv
import create_model
import corr_try

# Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument("-date", required=True)    # yyyy-mm
args = parser.parse_args()

date = args.date



# import models
models = genfromtxt('outputs/models/merged_2018-01.csv', delimiter=',')
models = np.transpose(models)

# Dates
dates_tHour = pd.date_range(
    pd.to_datetime(date),
    pd.to_datetime(date) + pd.DateOffset(months=6),
    closed='left',
    freq='1H'
    )

# observations
obs_H2, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO")   # could add st dev

# import correlation
#corr = corr_try.correlation(date)

# try without function

obs = obs_H2
        
    #remove nans
for i in range(0,len(obs)):
    if np.isnan(obs[i]):
        obs[i]=523.024
        
obs_H2, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO")   # could add st dev

    
    # create corr list
corr = np.ndarray(shape = (25))
i = 0
while i < 25:
    corr[i], _ = pearsonr(obs, models[i])
    i += 1



# Plots

fig = {}
ax = {}
color = cm.viridis(np.linspace(0, 1, 25))

fig_param = {
    'w': 6, 'h': 3,
    'px0': 0.80, 'py0': 0.50,
    'pw': 5.15, 'ph': 2.45,
    'ylblx': 0.05, 'ylbly': 1.5,  # left, centre aligned
    'fontsize1': 4,
    'fontsize' : 8,
    }
plt.close('all')



# Concentration
fig['main'] = plt.figure(figsize=(fig_param['w'], fig_param['h']), dpi=300)
for i in ['H2']:
    fig['main'].clf()
    ax['main'] =  name_qch4_couple.plot_h2.generic(
        fig=fig['main'],
        idata={

            'model1': [
                'line',
                [dates_tHour, np.array(models[0]), '-'],
                {'c': color[0], 'lw': 0.5, 'label': 'M1: {:.2f}'.format(corr[0])}
                ],
            'model2': [
                'line',
                [dates_tHour, np.array(models[1]), '-'],
                {'c': color[1], 'lw': 0.5, 'label': 'M2: {:.2f}'.format(corr[1])}
                ],
            'model3': [
                'line',
                [dates_tHour, np.array(models[2]), '-'],
                {'c': color[2], 'lw': 0.5, 'label': 'M3: {:.2f}'.format(corr[2])}
                ],
            'model4': [
                'line',
                [dates_tHour, np.array(models[3]), '-'],
                {'c': color[3], 'lw': 0.5, 'label': 'M4: {:.2f}'.format(corr[3])}
                ],
            'model5': [
                'line',
                [dates_tHour, np.array(models[4]), '-'],
                {'c': color[4], 'lw': 0.5, 'label': 'M5: {:.2f}'.format(corr[4])}
                ],
            'model6': [
                'line',
                [dates_tHour, np.array(models[5]), '-'],
                {'c': color[5], 'lw': 0.5, 'label': 'M6: {:.2f}'.format(corr[5])}
                ],
            'model7': [
                'line',
                [dates_tHour, np.array(models[6]), '-'],
                {'c': color[6], 'lw': 0.5, 'label': 'M7: {:.2f}'.format(corr[6])}
                ],
            'model8': [
                'line',
                [dates_tHour, np.array(models[7]), '-'],
                {'c': color[7], 'lw': 0.5, 'label': 'M8: {:.2f}'.format(corr[7])}
                ],
            'model9': [
                'line',
                [dates_tHour, np.array(models[8]), '-'],
                {'c': color[8], 'lw': 0.5, 'label': 'M9: {:.2f}'.format(corr[8])}
                ],
            'model10': [
                'line',
                [dates_tHour, np.array(models[9]), '-'],
                {'c': color[9], 'lw': 0.5, 'label': 'M10: {:.2f}'.format(corr[9])}
                ],
            'model11': [
                'line',
                [dates_tHour, np.array(models[10]), '-'],
                {'c': color[10], 'lw': 0.5, 'label': 'M11: {:.2f}'.format(corr[10])}
                ],
            'model12': [
                'line',
                [dates_tHour, np.array(models[11]), '-'],
                {'c': color[11], 'lw': 0.5, 'label': 'M12: {:.2f}'.format(corr[11])}
                ],
            'model13': [
                'line',
                [dates_tHour, np.array(models[12]), '-'],
                {'c': color[12], 'lw': 0.5, 'label': 'M13: {:.2f}'.format(corr[12])}
                ],
            'model14': [
                'line',
                [dates_tHour, np.array(models[13]), '-'],
                {'c': color[13], 'lw': 0.5, 'label': 'M14: {:.2f}'.format(corr[13])}
                ],
            'model15': [
                'line',
                [dates_tHour, np.array(models[14]), '-'],
                {'c': color[14], 'lw': 0.5, 'label': 'M15: {:.2f}'.format(corr[14])}
                ],
            'model16': [
                'line',
                [dates_tHour, np.array(models[15]), '-'],
                {'c': color[15], 'lw': 0.5, 'label': 'M16: {:.2f}'.format(corr[15])}
                ],
            'model17': [
                'line',
                [dates_tHour, np.array(models[16]), '-'],
                {'c': color[16], 'lw': 0.5, 'label': 'M17: {:.2f}'.format(corr[16])}
                ],
            'model18': [
                'line',
                [dates_tHour, np.array(models[17]), '-'],
                {'c': color[17], 'lw': 0.5, 'label': 'M18: {:.2f}'.format(corr[17])}
                ],
            'model19': [
                'line',
                [dates_tHour, np.array(models[18]), '-'],
                {'c': color[18], 'lw': 0.5, 'label': 'M19: {:.2f}'.format(corr[18])}
                ],
            'model20': [
                'line',
                [dates_tHour, np.array(models[19]), '-'],
                {'c': color[19], 'lw': 0.5, 'label': 'M20: {:.2f}'.format(corr[19])}
                ],
            'model21': [
                'line',
                [dates_tHour, np.array(models[20]), '-'],
                {'c': color[20], 'lw': 0.5, 'label': 'M21: {:.2f}'.format(corr[20])}
                ],
            'model22': [
                'line',
                [dates_tHour, np.array(models[21]), '-'],
                {'c': color[21], 'lw': 0.5, 'label': 'M22: {:.2f}'.format(corr[21])}
                ],
            'model23': [
                'line',
                [dates_tHour, np.array(models[22]), '-'],
                {'c': color[22], 'lw': 0.5, 'label': 'M23: {:.2f}'.format(corr[22])}
                ],
            'model24': [
                'line',
                [dates_tHour, np.array(models[23]), '-'],
                {'c': color[23], 'lw': 0.5, 'label': 'M24: {:.2f}'.format(corr[23])}
                ],
            'model25': [
                'line',
                [dates_tHour, np.array(models[24]), '-'],
                {'c': color[24], 'lw': 0.5, 'label': 'M25: {:.2f}'.format(corr[24])}
                ],
            'wao': [
                'line',
                [dates_tHour, np.array(obs_H2), '-'],
                {'c': '#FF0000', 'lw': 0.5, 'label': 'Observed WAO'}
                ],
            
              },

        texts=[
            {
                'x': fig_param['ylblx'] / fig_param['w'],
                'y': fig_param['ylbly'] / fig_param['h'],
                's': (
                    u'$\chi$ H$_{2}$ (nmol mol$^{-1}$)'
                    ),
                'ha': 'left', 'va': 'center',
                'size': fig_param['fontsize'], 'rotation': 90
                }
            ],
        xlim=[
            pd.to_datetime(date), 
            pd.to_datetime(date) + pd.DateOffset(months=6), 
            ],
        ylim=(
            [400., 660.]
            ),
        yticks=(
            np.arange(400., 660., 40.)
            ),
        tick_fontsize=fig_param['fontsize'],
        loc_plot=[
            fig_param['px0'] / fig_param['w'],
            fig_param['py0'] / fig_param['h'],
            fig_param['pw'] / fig_param['w'],
            fig_param['ph'] / fig_param['h']
            ],
        xtick_params=[
            True,
            mdates.DateFormatter('%Y-%m-%d'),
            mdates.MonthLocator(bymonthday=1),
            ]
        )
    for l in ax['main'].get_xticklabels():
        l.set_ha("right")
        l.set_rotation(30)
    ax['main'].legend(
        loc='upper right', ncol=7, fontsize=fig_param['fontsize1']
        )
    fig['main'].savefig(f'outputs/models/models_{date}_period_{i}.png')
