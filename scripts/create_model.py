# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 04:47:31 2022

@author: filot
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





def run_model(date):
    factor = 0
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
    
                  
    
    # Dilution matrix - H2 wao
    Dfile_H2 = (
        'inputs/footprints_wao/'
        f'WAO-20magl_UKV_EUROPE_{date_nodash}.nc'
        )
    with xr.open_dataset(Dfile_H2) as ds_read:
        with ds_read.load() as Din:
            D_wao = Din.fp.transpose('time', 'lat', 'lon').values
    
    
    # baseline
    if factor == 1:
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
                chi0 = obs_H2 + (D_wao * Q_b).sum((1, 2)) / M_H2 * 1e9
            mod = pd.Series(
                chi0 + (D_wao * Q).sum((1, 2)) / M_H2 * 1e9,
                index=dates_tHour
                )
            models[nom,] = mod.values
            nom += 1
    modelT = np.transpose(models)
    np.savetxt(f'outputs/models/models_changedbaseline_{date}.csv', modelT, delimiter=",")
  #  return models


