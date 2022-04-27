# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 04:47:31 2022

Creating the tweaked modelled 'baselines'

@author: Filoteea Moldovan
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







# =============================================================================
# function for creating tweaked modelled 'baselines' from Section 2.9.
# =============================================================================

def run_model_bas(date, factor):
    date_nodash = date.replace('-', '')    
    # Dates
    dates_tHour = pd.date_range(
        pd.to_datetime(date),
        pd.to_datetime(date) + pd.DateOffset(months=1),
        closed='left',
        freq='1H'
        )
    
    M_H2 = 2.016
                 
    
    if factor == 1:  #model basline for MHD
        Dfile_H2 = (
            'inputs/baseline/footprints_mhd/'
            f'MHD-10magl_UKV_EUROPE_{date_nodash}.nc'
            )
        with xr.open_dataset(Dfile_H2) as ds_read:
            with ds_read.load() as Din:
                D = Din.fp.transpose('time', 'lat', 'lon').values
        chi, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "MHD_10magl", 0)  # read observations 
        for i in range(0,len(chi)): #interpolate missing datapoints
            if np.isnan(chi[i]):
                chi[i] = 523.024   # data average

    else:   #model baseline for WAO
        Dfile_H2 = (
            'inputs/footprints_wao/'
            f'WAO-20magl_UKV_EUROPE_{date_nodash}.nc'
            )
        with xr.open_dataset(Dfile_H2) as ds_read:
            with ds_read.load() as Din:
                D = Din.fp.transpose('time', 'lat', 'lon').values
        # read WAO observations
        chi, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO", 0)  #read observations 
        for i in range(0,len(chi)): #interpolate missing datapoints
            if np.isnan(chi[i]):
                chi[i] = 510.400  #data average        
    
    # create all tweaked baselines simultaneously 
    '''
    sink = [0.01, 0.1, 1, 2.5, 5]
    emission_coeff = [2, 1.5, 1, 0.7, 0.5]
    models = np.ndarray(shape = (25, len(dates_tHour)))
    nom = 0
    for i in emission_coeff:
        for j in sink:
            Q = chem_co.read_Qsink(dates_tHour, i, j, 0)                   
            mod = pd.Series(
                chi - (D * Q).sum((1, 2)) / M_H2 * 1e9,
                index=dates_tHour
                )
            models[nom,] = mod.values
            nom += 1
    modelT = np.transpose(models)
    np.savetxt(f'outputs/models/baselines/higher_emm/baseline_mhd_2xemm_{date}.csv', modelT, delimiter=",")
    '''
    
    # create final model 'baseline'
    '''
    Q = chem_co.read_Qsink(dates_tHour, 1, 1, 0)                   
    mod = pd.Series(
        chi - (D * Q).sum((1, 2)) / M_H2 * 1e9,
        index=dates_tHour
        )   
           
    pd.concat([
    pd.Series(mod, index=dates_tHour, name='h2_ppb'),
    ], axis=1).to_csv(f'outputs/scenarios/bas_wao_{date}.csv') 
  #  return models
  '''

# create tweaked modesl for each month

dates = ['2018-01', '2018-02', '2018-03', '2018-04', '2018-05', '2018-06',
         '2018-07', '2018-08', '2018-09', '2018-10', '2018-11', '2018-12']

x = 1 # for MHD
x = 2 # for WAO

for i in dates:
    run_model_bas(i, x)





# function for creating 25 tweaked models (not included in dissertation)

def run_model(date):
    factor = 2
    date_nodash = date.replace('-', '')    
    # Dates
    dates_tHour = pd.date_range(
        pd.to_datetime(date),
        pd.to_datetime(date) + pd.DateOffset(months=1),
        closed='left',
        freq='1H'
        )
    
    M_H2 = 2.016
    # Dilution matrix - H2 MHD
    Dfile_H2 = (
        'inputs/baseline/footprints_mhd/'
        f'MHD-10magl_UKV_EUROPE_{date_nodash}.nc'
        )
    with xr.open_dataset(Dfile_H2) as ds_read:
        with ds_read.load() as Din:
            D_mhd = Din.fp.transpose('time', 'lat', 'lon').values
                 
    
    # Dilution matrix - H2 wao
    Dfile_H2 = (
        'inputs/footprints_wao/'
        f'WAO-20magl_UKV_EUROPE_{date_nodash}.nc'
        )
    with xr.open_dataset(Dfile_H2) as ds_read:
        with ds_read.load() as Din:
            D_wao = Din.fp.transpose('time', 'lat', 'lon').values
    
    
    # baseline
    if factor == 2:
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
        
        bas = read_baseline(dates_tHour)
   # elif factor == 2:
        chi0, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "MHD_10magl")   # could add st dev

        for i in range(0,len(chi0)):
            if np.isnan(chi0[i]):
                chi0[i] = bas[i]     #523.024
    else:
        # read WAO observations
        obs_H2, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO")   # could add st dev

        #remove nans
        for i in range(0,len(obs_H2)):
            if np.isnan(obs_H2[i]):
                obs_H2[i]=523.024

      
    
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
            if factor == 0:
                Q_b = chem_co.read_Qsink(dates_tHour, 0, j, 0)
                chi0_proc = chi0 - (D_mhd * Q_b).sum((1, 2)) / M_H2 * 1e9
                           # save new baseline
                if i == 2 and j == 1:
                    np.savetxt(f'outputs/models/mhd_bas/{date}_chi0p.csv', chi0_proc, delimiter=",")
                    
            else:
                chi0_proc = chi0 - (D_mhd * Q).sum((1, 2)) / M_H2 * 1e9
                if i == 1.5 and j == 1:
                    np.savetxt(f'outputs/models/mhd_bas/{date}_chi0p_obs.csv', chi0_proc, delimiter=",")

                                   
            mod = pd.Series(
                chi0_proc + (D_wao * Q).sum((1, 2)) / M_H2 * 1e9,
                index=dates_tHour
                )
            models[nom,] = mod.values
            nom += 1
    modelT = np.transpose(models)
    np.savetxt(f'outputs/models/mhd_bas/models_mhd_obs_{date}.csv', modelT, delimiter=",")