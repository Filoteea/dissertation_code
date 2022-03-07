
"""
Created on Wed Jan 26 13:47:47 2022

@author: filot
Adapted from Edward Chung
"""

r"""

H2 Forward Model verification

"""
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


# Third party imports
from collections import OrderedDict
from datetime import datetime

# Semi-local imports
import name_qch4_couple.io
import name_qch4_couple.name
import name_qch4_couple.plot_h2


# Local imports
import routines
import chem_co


# Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument("-date", required=True)    # yyyy-mm
parser.add_argument("-odir", required=True)
args = parser.parse_args()

date = args.date
odir = args.odir

date_nodash = date.replace('-', '')

# Dates
dates_tHour = pd.date_range(
    pd.to_datetime(date),
    pd.to_datetime(date) + pd.DateOffset(months=1),
    closed='left',
    freq='1H'
    )

# Grid
grid_info = routines.define_grid()
inv_reg_map0 = grid_info['inv_reg_map']
nlat = grid_info['nlat']
nlon = grid_info['nlon']
area = grid_info['area']
grid_centre = grid_info['grid_centre']
grid_vertex = grid_info['grid_vertex']
inv_reg_uniq = grid_info['inv_reg_uniq']

M_H2 = 2.016


# Emission scenarios
#Q = chem_co.read_Q(dates_tHour)
              
# Scenario emission Qs

'''
conversion documentation:
    1 kg CH4 to 1 kg H2: 0.368
    H2 demand (by 2018 natural gas consumption in 2018):
        2035 -> 0.147
        2050 -> 0.687
        ideal -> 1
    leakage: 
        1% -> 0.01
        10% -> 0.1
        50% -> 0.5
'''
#factor = 0.368
#Qs = chem_co.read_Qs(dates_tHour, factor)

# Dilution matrix - H2 wao
Dfile_H2 = (
    'inputs/footprints_wao/'
    f'WAO-20magl_UKV_EUROPE_{date_nodash}.nc'
    )
with xr.open_dataset(Dfile_H2) as ds_read:
    with ds_read.load() as Din:
        D_wao = Din.fp.transpose('time', 'lat', 'lon').values

'''
# Dilution matrix - H2 mhd
Dfile_H2 = (
    'inputs/baseline/footprints_mhd/'
    f'MHD-10magl_UKV_EUROPE_{date_nodash}.nc'
    )
with xr.open_dataset(Dfile_H2) as ds_read:
    with ds_read.load() as Din:
        D_mhd = Din.fp.transpose('time', 'lat', 'lon').values
'''

# baseline

def read_baseline(timestamps):
    date = timestamps[0].strftime('%Y-%m')
    year = timestamps[0].strftime('%Y')
    chi0file = (
                'outputs/baseline/baseline-MHD_10magl-h2-2018.nc'
                )
    with xr.open_dataset(chi0file) as ds_read:     #put as
        with ds_read.load() as ds:
                chi0 = ds.chi_H2.sel(time=date).to_series()
                
    return chi0

chi0 = read_baseline(dates_tHour)
# average+bas = 523.024

# models

'''
# wao
mod_H2_bas = pd.Series(
    chi0 + (D_wao * Q).sum((1, 2)) / M_H2 * 1e9,
    index=dates_tHour
    )
'''

'''
# mhd 
mod_mhd = pd.Series(
    chi0 + (D_mhd * Q).sum((1, 2)) / M_H2 * 1e9,
    index=dates_tHour
    )
'''    

# create 25 models
'''
# mod scenarios
  25 models
  - 5 sink scenarios:
      -- 0.01
      -- 0.1
      -- 1
      -- 10
      -- 50
      
  - 5 current emission coeff:
      -- 2
      -- 1.5
      -- 1
      -- 0.7
      -- 0.5
'''
sink = [0.01, 0.1, 1, 2.5, 5]
emission_coeff = [2, 1.5, 1, 0.7, 0.5]
models = np.ndarray(shape = (25, len(dates_tHour)))
nom = 0

for i in emission_coeff:
    for j in sink:
        Q = chem_co.read_Qsink(dates_tHour, i, j, 0)
        mod = pd.Series(
            chi0 + (D_wao * Q).sum((1, 2)) / M_H2 * 1e9,
            index=dates_tHour
            )
        models[nom,] = mod.values
        nom += 1

'''
>>> numpy.zeros((3, 5))
    array([[ 0.,  0.,  0.,  0.,  0.],
   [ 0.,  0.,  0.,  0.,  0.],
   [ 0.,  0.,  0.,  0.,  0.]])
'''



# observations
obs_H2, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO")   # could add st dev

# obs_H2_mhd, sigma_obs_H2_mhd = chem_co.read_obs(dates_tHour, "MHD_10magl")   # could add st dev


# calculate correlation
#corr, _ = spearmanr(models[0], obs_H2)



# to save 

modelT = np.transpose(models)
np.savetxt(os.path.join(odir, f'models_{date}.csv'), modelT, delimiter=",")

# Plots

fig = {}
ax = {}
color = cm.viridis(np.linspace(0, 1, 25))
colours = {
    'obs': '#000000',
    'mo1': '#8888FF',
    'bas': '#886600',
    'mo2': "#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)]),
    }
fig_param = {
    'w': 6, 'h': 3,
    'px0': 0.80, 'py0': 0.50,
    'pw': 5.15, 'ph': 2.45,
    'ylblx': 0.05, 'ylbly': 1.5,  # left, centre aligned
    'fontsize': 6,
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
                {'c': color[0], 'lw': 0.5, 'label': 'M1'}
                ],
            'model2': [
                'line',
                [dates_tHour, np.array(models[1]), '-'],
                {'c': color[1], 'lw': 0.5, 'label': 'M2'}
                ],
            'model3': [
                'line',
                [dates_tHour, np.array(models[2]), '-'],
                {'c': color[2], 'lw': 0.5, 'label': 'M3'}
                ],
            'model4': [
                'line',
                [dates_tHour, np.array(models[3]), '-'],
                {'c': color[3], 'lw': 0.5, 'label': 'M4'}
                ],
            'model5': [
                'line',
                [dates_tHour, np.array(models[4]), '-'],
                {'c': color[4], 'lw': 0.5, 'label': 'M5'}
                ],
            'model6': [
                'line',
                [dates_tHour, np.array(models[5]), '-'],
                {'c': color[5], 'lw': 0.5, 'label': 'M6'}
                ],
            'model7': [
                'line',
                [dates_tHour, np.array(models[6]), '-'],
                {'c': color[6], 'lw': 0.5, 'label': 'M7'}
                ],
            'model8': [
                'line',
                [dates_tHour, np.array(models[7]), '-'],
                {'c': color[7], 'lw': 0.5, 'label': 'M8'}
                ],
            'model9': [
                'line',
                [dates_tHour, np.array(models[8]), '-'],
                {'c': color[8], 'lw': 0.5, 'label': 'M9'}
                ],
            'model10': [
                'line',
                [dates_tHour, np.array(models[9]), '-'],
                {'c': color[9], 'lw': 0.5, 'label': 'M10'}
                ],
            'model11': [
                'line',
                [dates_tHour, np.array(models[10]), '-'],
                {'c': color[10], 'lw': 0.5, 'label': 'M11'}
                ],
            'model12': [
                'line',
                [dates_tHour, np.array(models[11]), '-'],
                {'c': color[11], 'lw': 0.5, 'label': 'M12'}
                ],
            'model13': [
                'line',
                [dates_tHour, np.array(models[12]), '-'],
                {'c': color[12], 'lw': 0.5, 'label': 'M13'}
                ],
            'model14': [
                'line',
                [dates_tHour, np.array(models[13]), '-'],
                {'c': color[13], 'lw': 0.5, 'label': 'M14'}
                ],
            'model15': [
                'line',
                [dates_tHour, np.array(models[14]), '-'],
                {'c': color[14], 'lw': 0.5, 'label': 'M15'}
                ],
            'model16': [
                'line',
                [dates_tHour, np.array(models[15]), '-'],
                {'c': color[15], 'lw': 0.5, 'label': 'M16'}
                ],
            'model17': [
                'line',
                [dates_tHour, np.array(models[16]), '-'],
                {'c': color[16], 'lw': 0.5, 'label': 'M17'}
                ],
            'model18': [
                'line',
                [dates_tHour, np.array(models[17]), '-'],
                {'c': color[17], 'lw': 0.5, 'label': 'M18'}
                ],
            'model19': [
                'line',
                [dates_tHour, np.array(models[18]), '-'],
                {'c': color[18], 'lw': 0.5, 'label': 'M19'}
                ],
            'model20': [
                'line',
                [dates_tHour, np.array(models[19]), '-'],
                {'c': color[19], 'lw': 0.5, 'label': 'M20'}
                ],
            'model21': [
                'line',
                [dates_tHour, np.array(models[20]), '-'],
                {'c': color[20], 'lw': 0.5, 'label': 'M21'}
                ],
            'model22': [
                'line',
                [dates_tHour, np.array(models[21]), '-'],
                {'c': color[21], 'lw': 0.5, 'label': 'M22'}
                ],
            'model23': [
                'line',
                [dates_tHour, np.array(models[22]), '-'],
                {'c': color[22], 'lw': 0.5, 'label': 'M23'}
                ],
            'model24': [
                'line',
                [dates_tHour, np.array(models[23]), '-'],
                {'c': color[23], 'lw': 0.5, 'label': 'M24'}
                ],
            'model25': [
                'line',
                [dates_tHour, np.array(models[24]), '-'],
                {'c': color[24], 'lw': 0.5, 'label': 'M25'}
                ],
            'wao': [
                'line',
                [dates_tHour, np.array(obs_H2), '-'],
                {'c': '#FF0000', 'lw': 0.5, 'label': 'Observed WAO'}
                ],

            #'mhd': [
            #    'line',
            #    [obs_H2_mhd.index, np.array(obs_H2_mhd), '-'],
            #    {'c': '#0139FF', 'lw': 0.5, 'label': 'Measured MHD'}
            #    ],
            #obs2': [
            #    'line',
            #    [mod_H2_bas.index, np.array(mod_H2_bas), '-'],
            #    {'c': colours['bas'], 'lw': 0.5, 'label': 'Modelled with baseline'}
            #    ],     
            #'bas': [
            #    'line',
             #   [chi0.index, chi0.values, '-'],
             #   {'c': colours['bas'], 'lw': 0.5, 'label': 'Baseline'}
             #   ],              
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
            pd.to_datetime(date) + pd.DateOffset(months=1), 
            ],
        ylim=(
            [430., 650.]
            ),
        yticks=(
            np.arange(430., 650., 20.)
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
            mdates.WeekdayLocator(byweekday=6),
            ]
        )
    for l in ax['main'].get_xticklabels():
        l.set_ha("right")
        l.set_rotation(30)
    ax['main'].legend(
        loc='upper right', ncol=7, fontsize=fig_param['fontsize']
        )
    fig['main'].savefig(os.path.join(odir, f'model_trial_{date}_{i}.png'))
