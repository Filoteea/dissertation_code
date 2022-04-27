# -*- coding: utf-8 -*-
"""
Created on Thur Mar 24 18:01:48 2022

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
from scipy.stats import pearsonr
from numpy import genfromtxt


# Third party imports
from collections import OrderedDict
from datetime import datetime

# Semi-local imports
import name_qch4_couple.io
import name_qch4_couple.name



# Local imports
import chem_co


# Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument("-date", required=True)    # yyyy-mm
args = parser.parse_args()

date = args.date



# function that runs the model
def run(date, factor):
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
        
# =============================================================================
# create emission scenarios for HFD
# =============================================================================

        # Dilution matrix - H2 wao
    Dfile_H2 = (
        'inputs/footprints_wao/'
        f'WAO-20magl_UKV_EUROPE_{date_nodash}.nc'
        )
    with xr.open_dataset(Dfile_H2) as ds_read:
        with ds_read.load() as Din:
            D_wao = Din.fp.transpose('time', 'lat', 'lon').values
    
    
# =============================================================================
# create emission scenarios for HFD
# =============================================================================
    
        # Dilution matrix - H2 HFD
   # Dfile_H2 = (
   #     'inputs/footprints_hfd/'
   #     f'HFD-100magl_UKV_EUROPE_{date_nodash}.nc'
   #     )
   # with xr.open_dataset(Dfile_H2) as ds_read:
   #     with ds_read.load() as Din:
   #         D_wao = Din.fp.transpose('time', 'lat', 'lon').values
        
        
        # calculate modelled 'baseline' (chi0)
    chi1, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "MHD_10magl", 0)   # could add st dev
    for i in range(0,len(chi1)):
        if np.isnan(chi1[i]):
            chi1[i] = 523.024   # data average
    Q_0 = chem_co.read_Qsink(dates_tHour, 1, 1, 0) 
                        
    chi0 = pd.Series(
        chi1 - (D_mhd * Q_0).sum((1, 2)) / M_H2 * 1e9,
        index = dates_tHour
        )
    
       # calculate projected H2 concentrations (factor = Emission Coefficient; Table 2.4)        
    Q = chem_co.read_Qsink(dates_tHour, 1, 1, factor)
    mod = pd.Series(
        chi0 + (D_wao * Q).sum((1, 2)) / M_H2 * 1e9,
        index = dates_tHour
       )
    
    pd.concat([
        pd.Series(mod, index=dates_tHour, name='h2_ppb'),
        ], axis=1).to_csv(f'outputs/scenarios/new_scenario/scenario_wao_{factor}_{date}.csv')



#run model for all months and emission factors
dates = ['2018-01', '2018-02', '2018-03', '2018-04', '2018-05', '2018-06',
         '2018-07', '2018-08', '2018-09', '2018-10', '2018-11', '2018-12']

factor = [0, 0.035, 0.165, 0.241]

for j in factor:
    for i in dates:
        run(i, j)



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
        1000% -> 1
        
    final values:
        2035: 0.0005, 0.0054, 0.027, 0.054
        2050: 0.0025, 0.0252, 0.126, 0.252
        ideal:0.0036, 0.0368, 0.184, 0.368
        
    new final: 
        2018: 0
        2035: 0.035
        2050: 0.165
        ideal: 0.241
'''
