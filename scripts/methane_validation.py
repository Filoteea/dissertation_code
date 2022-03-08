r"""

Radon/CH4 Forward Model

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

# Third party imports
from collections import OrderedDict
from datetime import datetime
from sklearn import linear_model

# Semi-local imports
import name_qch4_couple.io
import name_qch4_couple.name
import name_qch4_couple.plot_h2

# Local imports
import routines
import chem_ch4_validation
import chem_co


# Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument("-date", required=True)
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

# Standard atmospheric conditions
p_std = 1013.25
T_std = 15

# =============================================================================
#   CH4
# =============================================================================
M_CH4 = 16.043  # g mol-1 - IUPAC
M_H2 = 2.016

# Q - CH4
Q = chem_ch4_validation.read_Q(dates_tHour)


# Dilution matrix - CH4 mhd
Dfile_H2 = (
    'inputs/footprints_hfd/'
    f'HFD-100magl_UKV_EUROPE_{date_nodash}.nc'
    )
with xr.open_dataset(Dfile_H2) as ds_read:
    with ds_read.load() as Din:
        D = Din.fp.transpose('time', 'lat', 'lon').values

# baseline
def read_baseline(timestamps):
    date = timestamps[0].strftime('%Y-%m')
    year = timestamps[0].strftime('%Y')
    chi0file = (
                'outputs/validation_ch4/baseline-MHD_10magl-ch4-2018.nc'
                )
    with xr.open_dataset(chi0file) as ds_read:     #put as
        with ds_read.load() as ds:
                chi0 = ds.chi_CH4.sel(time=date).to_series()
                
    return chi0

chi0 = read_baseline(dates_tHour)

'''
# modelled methane mhd
mod_CH4 = pd.Series(
    chi0.values + (D * Q).sum((1, 2)) / M_CH4 * 1e9,
    index=dates_tHour
    )
'''
# for conversion
mod_CH4 = pd.Series(
    (D * Q).sum((1, 2)) / M_CH4 * 1e9,
    index=dates_tHour
    )
mod_H2 = pd.Series(
    (D * Q).sum((1, 2)) / M_H2 * 1e9,
    index=dates_tHour
    )
# obs mhd
#obs_ch4_mhd, sigma_obs_H2_mhd = chem_ch4_validation.read_obs(dates_tHour, "MHD_10magl")   # could add st dev

# to save 
#pd.concat([
#    pd.Series(mod_CH4, index=dates_tHour, name='chi_CH4'),
#    ], axis=1).to_csv(os.path.join(odir, f'hfd_ch4_{date}.csv'))

#'''
# Plots
fig = {}
ax = {}
colours = {
    'obs_mhd': '#000000',
    'mod_ch4': '#0000FF',
    }
fig_param = {
    'w': 6, 'h': 3,
    'px0': 0.80, 'py0': 0.50,
    'pw': 5.15, 'ph': 2.45,
    'ylblx': 0.05, 'ylbly': 1.5,  # left, centre aligned
    'fontsize': 8,
    }
plt.close('all')

# Concentration
fig['main'] = plt.figure(figsize=(fig_param['w'], fig_param['h']), dpi=300)
for i in ['CH4']:
    fig['main'].clf()
    ax['main'] =  name_qch4_couple.plot_h2.generic(
        fig=fig['main'],
        idata={
            #'hfd': [
            #    'line',
            #    [obs_ch4_mhd.index, np.array(obs_ch4_mhd), '-'],
            #    {'c': colours['obs_mhd'], 'lw': 0.5, 'label': 'Measured HFD'}
            #    ],
            'mod_h2': [
                'line',
                [mod_H2.index, np.array(mod_H2), '-'],
                {'c': colours['obs_mhd'], 'lw': 0.5, 'label': u'Converted to H$_{2}$'}
                ],
            'mod_ch4': [
                'line',
                [mod_CH4.index, np.array(mod_CH4), '-'],
                {'c': colours['mod_ch4'], 'lw': 0.5, 'label': u'Modelled CH$_{4}$'}
                ],
            },

        texts=[
            {
                'x': fig_param['ylblx'] / fig_param['w'],
                'y': fig_param['ylbly'] / fig_param['h'],
                's': (
                    u'Concentration (nmol mol$^{-1}$)'
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
            [0., 2000.]
            ),
        yticks=(
            np.arange(0., 2000., 200.)
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
        loc='upper right', ncol=4, fontsize=fig_param['fontsize']
        )
    fig['main'].savefig(os.path.join(odir, f'conversion_{i}.png'))
#'''
