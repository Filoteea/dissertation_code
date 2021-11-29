r"""

Baseline Calculation

"""
# Standard Library imports
import argparse
import datetime
import h5py
import json
import matplotlib.colors
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import os
import pandas as pd
import re
import scipy.optimize
import warnings
import sys
import xarray as xr

# Third party imports
from collections import OrderedDict

# Semi-local imports
import name_qch4_couple.io
import name_qch4_couple.name
import name_qch4_couple.plot
import name_qch4_couple.region_EU
import name_qch4_couple.routines
import name_qch4_couple.util

# Local imports
import routines

# =============================================================================
#   Settings
# =============================================================================
# Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument("-site", required=True)
parser.add_argument("-year", required=True, type=int)
parser.add_argument("-window1", required=True, type=int)
parser.add_argument("-window2", required=True, type=int)
parser.add_argument("-force", required=True, type=int)
parser.add_argument("-odir", required=True)
args = parser.parse_args()

site = args.site
year = args.year
window1 = args.window1
window2 = args.window2
force_compute = bool(args.force)
odir = args.odir

mode = 'all'
site_ref = 'mhd'  # HDF reference key
p = 101325  # Pa
T = 288.15  # K

ifile_fil = f'filtered_result_{mode}_{year}.hdf'
ofile_fil = f'filtered_result_{mode}_{year}.hdf'
ifile_con = f'filtered_condition_{mode}_{year}.hdf'
ofile_con = f'filtered_condition_{mode}_{year}.hdf'

ofile_fin = f'{site}_{mode}_{year}.hdf'
ofile_fin_nc = f'{site}_{mode}_{year}.nc'

field_dirs = OrderedDict()
particle_dirs = OrderedDict()
obs_files = OrderedDict()
long_names = OrderedDict()
locations = OrderedDict()
field_prefix = OrderedDict()
particle_prefix = OrderedDict()
with open(f'/home/ec5/name_qch4_couple/stations/{site}.json', 'r') as f:
    st_info = json.load(f)
    field_dirs[site] = OrderedDict(st_info['mix_dir'])
    particle_dirs[site] = OrderedDict(st_info['particle_dir'])
    obs_files[site] = OrderedDict(st_info['obs_files'])
    long_names[site] = st_info['long_name']
    locations[site] = st_info['location']
    field_prefix[site] = OrderedDict(st_info['mix_prefix'])
    particle_prefix[site] = OrderedDict(st_info['particle_prefix'])

# =============================================================================
#   Pre-processing
# =============================================================================
print(f'Initialising')
print(f'    site = {site}')
print(f'    year = {year}')

# Grid
grid_info = routines.define_grid()
dlat = grid_info['dlat']
nlat = grid_info['nlat']
nlon = grid_info['nlon']
area = grid_info['area']
grid_centre = grid_info['grid_centre']
grid_vertex = grid_info['grid_vertex']

# Dates
dt = pd.to_timedelta(window1//2, 'H')
dates_tHour = pd.date_range(
    pd.to_datetime(f'{year}') - pd.to_timedelta(window1//2, 'H'),
    pd.to_datetime(f'{int(year) + 1}') + pd.to_timedelta(window1//2, 'H'),
    freq='1H',
    closed='left'
    )
dates_tDay = dates_tHour.floor('D').unique()
hdf_key = (
    f'{site}_{year}'
    )

# Observations
obs, obs_err = routines.read_obs(
    sites={site: {k: 0 for k in ['chi', 'd_13C', 'd_2H']}},
    obs_files=obs_files,
    date_range=[dates_tDay[i].strftime('%Y-%m-%d') for i in [0, -1]],
    dates_tHour=dates_tHour
    )

# Filters
filter1 = {
    k: np.zeros((nlat, nlon)) for k in ['f_local', 'f_europe', 'f_south']
    }

# Local
radius = 2
st_idx = [np.abs(grid_centre[i] - l).argmin() for i, l in enumerate(locations[site])]
filter1['f_local'][
    st_idx[1]-radius:st_idx[1]+radius+1,
    st_idx[0]-radius:st_idx[0]+radius+1
    ] = 1.

# Europe
filter1['f_europe'][
    np.load('/home/ec5/name_qch4_couple/shape_files/europe.npy', 'r')
    ] = 1.

# Southern
lat_south = 30
filter1['f_south'][grid_centre[1] + 0.5*dlat <= lat_south] = 1.

# Population
# filter1['f_pop'] = 

# Threshold
filter_threshold = {k: 1. for k in filter1}
filter_threshold['f_local'] = 10
filter_threshold['f_europe'] = 5
filter_threshold['f_south'] = 2

# Dilution sensitivity limit
dsl = 3.4 * 1.e6 / p * 8.314 * T

# Domain borders
filter2 = {
    'd_bg': [1, 6, 7, 8],
    'd_bg2': [1, 6, 7, 8, 9, 14, 15, 16],  # default
    'd_bg3': [1, 6, 7, 8, 9, 14, 15, 16, 17]
    }
filter_threshold['d_bg'] = 0.8
filter_threshold['d_bg2'] = 0.8
filter_threshold['d_bg3'] = 0.8

# =============================================================================
#   Main Process
# =============================================================================
print("Processing Main Routine")
filters = list(filter1) + list(filter2)
flag_con = False

# Get conditions
if ifile_con and not force_compute and os.path.exists(os.path.join(odir, ifile_con)):
    flag_con = name_qch4_couple.util.hdf_key_check(
        os.path.join(odir, ifile_con), hdf_key
        )
if flag_con:
    print("    Using pre-existing filters")
    condition = pd.read_hdf(
        os.path.join(odir, ifile_con), hdf_key
        )
    idx_H = condition.index
    idx_d = idx_H.floor('D').unique()
else:
    print("    Applying filters")
    if mode == 'obs':
        idx_H = pd.concat(
                [v for k, v in obs[site].items()]
                ).sort_index().index.floor('H').unique()
    else:
        idx_H = dates_tHour 
    idx_d = idx_H.floor('D').unique()
    condition = pd.DataFrame(
            0,
            index=idx_H,
            columns=filters
            )
    # Fields files
    f_fields_paths = {
        site:
        pd.concat(
            [
                dates_tDay.strftime(
                    os.path.join(v[k], field_prefix[site][k])
                    ).to_series(index=dates_tDay).apply(np.vectorize(
                        lambda x: x if os.path.exists(x) else False
                        ))
                for k in sorted(v)
                ],
            axis=1
            ).replace(False, np.nan)
        for site, v in field_dirs.items()
        }
    f_fields_path = {
        site:
        v.bfill(axis=1).iloc[:, 0].replace(np.nan, False)
        for site, v in f_fields_paths.items()
        }
    # Particles files
    f_particles_paths = {
        site:
        pd.concat(
            [
                dates_tDay.strftime(
                    os.path.join(v[k], particle_prefix[site][k])
                    ).to_series(index=dates_tDay).apply(np.vectorize(
                        lambda x: x if os.path.exists(x) else np.nan
                        ))
                for k in sorted(v)
                ],
            axis=1
            ).replace(False, np.nan)
        for site, v in particle_dirs.items()
        }
    f_particles_path = {
        site:
        v.bfill(axis=1).iloc[:, 0].replace(np.nan, False)
        for site, v in f_particles_paths.items()
        }
    date_iter = idx_H.to_series().groupby(idx_H.date).apply(pd.DatetimeIndex)
    date_iter.index = pd.DatetimeIndex(date_iter.index)
    for date, hours in date_iter.iteritems():
        try:
            # Fields
            f_fields = f_fields_path[site].loc[date]
            dosage = name_qch4_couple.name.txt2array(
                name_qch4_couple.io.r_name(f_fields),
                [nlat, nlon], hours.hour, fill=0.0
                )
            # dilution = dosage * area / 3600.  # g m-3 s * m2 / g = m-1 s
            dilution = dosage * area / 3600.  # ppm s * m2 / g = ppm s m2 g-1
            # Particles
            f_particles = f_particles_path[site].loc[date]
            end_loc = name_qch4_couple.name.cat_end_loc(
                name_qch4_couple.io.r_name_particles(f_particles),
                grid_vertex, hours.hour, 3.e3, 8.e3
                )[0]
            end_loc_sel = end_loc[hours.hour]
            # Conditions
            condition.loc[hours, filter1] = pd.DataFrame(
                {
                    k: (dilution[hours.hour] * v).sum(axis=(1, 2))
                    for k, v in filter1.items()
                    },
                index=hours
                )
            condition.loc[hours, filter2] = pd.DataFrame(
                {
                    k: end_loc_sel[:, v].sum(1) / end_loc_sel.sum(1)
                    for k, v in filter2.items()
                    },
                index=hours
                )
        except:
            condition.loc[hours, filter1] = np.inf
            condition.loc[hours, filter2] = -np.inf
    condition.to_hdf(os.path.join(odir, ofile_con), hdf_key, mode='a')

# Check if baseline condition
filtered = pd.DataFrame(0, index=idx_H, columns=filters)
for k in filter1:
    threshold = filter_threshold[k]
    filtered.loc[condition[k] / dsl <= threshold, k] = 1.
for k in filter2:
    threshold = filter_threshold[k]
    filtered.loc[condition[k] >= threshold, k] = 1.
filtered.to_hdf(os.path.join(odir, ofile_fil), hdf_key, mode='a')

# Baseline
filters_fin = ['f_local', 'f_europe', 'f_south', 'd_bg2']
oidx = {k: v.index for k, v in obs[site].items()}
nidx = {k: v.floor('H') for k, v in oidx.items()}
obs_filter = {
    k: filtered.reindex(v).reset_index(drop=True)
    for k, v in nidx.items()
    }
for k, v in obs_filter.items():
    v.index = oidx[k]
baseline0 = {
    k: v[obs_filter[k][filters_fin].min(1) == 1.]
    for k, v in obs[site].items()
    }
baseline0_sigma = {
    k: v[obs_filter[k][filters_fin].min(1) == 1.]
    for k, v in obs_err[site].items()
    }

if not force_compute and os.path.exists(os.path.join(odir, ofile_fin_nc)):
    print("    Reading baselines")
    baseline2 = {}
    baseline2_var = {}
    with xr.open_dataset(os.path.join(odir, ofile_fin_nc)) as ds_read:
        with ds_read.load() as ds:
            for k in obs[site].keys():
                if not k in ds.variables:
                    continue
                baseline2[k] = ds[k].to_pandas()
            for k in obs[site].keys():
                if not f'var_{k}' in ds.variables:
                    continue
                baseline2_var[k] = ds[f'var_{k}'].to_pandas()
else:
    print("    Processing baselines")
    baseline1 = {k: pd.Series(np.nan, index=dates_tHour) for k in obs[site]}
    baseline1_var = {k: pd.Series(np.nan, index=dates_tHour) for k in obs[site]}
    warnings.simplefilter('ignore', np.RankWarning)
    check = pd.DataFrame(0., index=dates_tHour, columns=[])
    dt = pd.to_timedelta(f'{window1//2}H')
    for i in dates_tHour[window1//2+1:-(window1-window1//2+1)]:
        selected_data0 = {  # Select data within the default window
                k: v.loc[i - dt: i + dt]
                for k, v in baseline0.items()
                }
        nrow = {  # Minimum number of measurements before/after a given time
                k: 2 if i in v.index else 1
                #k: 3 if i in v.index else 2
                for k, v in selected_data0.items()
                }
        nbefore = {  # If there are at least the minimum number
                k: v[:i].size >= nrow[k]
                for k, v in selected_data0.items()
                }
        nafter = {  # If there are at least the minimum number
                k: v[i:].size >= nrow[k]
                for k, v in selected_data0.items()
                }
        date_range0 = {  # Extend window if necessary to include at least 1 point
                k:
                i - dt
                if nbefore[k] else
                baseline0[k].loc[:i].iloc[-nrow[k]:].first_valid_index()
                for k, v in selected_data0.items()
                }
        date_range1 = {  # Extend window if necessary to include at least 1 point
                k:
                i + dt
                if nafter[k] else
                baseline0[k].loc[i:].iloc[:nrow[k]].last_valid_index()
                for k, v in selected_data0.items()
                }
        selected_data = {  # Select data within the new window
                k: v.loc[date_range0[k]:date_range1[k]]
                for k, v in baseline0.items()
                }
        selected_data_sigma = {
                k: v.loc[date_range0[k]:date_range1[k]]
                for k, v in baseline0_sigma.items()
                }
        for k, v in selected_data.items():
            """
            3 measurement points: quadratic fit
            2 measurement points: linear fit
            1 measurement points: the value itself
            0 measurement points: NaN
            """
            if v.size >= 2 and None not in [date_range0[k], date_range1[k]]:
                order = min([v.size-1, 2]) if nbefore[k] and nafter[k] else 1
                td = (v.index - i).to_series() / pd.Timedelta(hours=1)
                td_i = 0  # (i - date_range0[k]) / pd.Timedelta(hours=1)
                fit = scipy.optimize.curve_fit(
                    (
                        (lambda x, a, b, c: a + b*x + c*x**2)
                        if order == 2 else
                        (lambda x, a, b: a + b*x)
                        ),
                    np.array(td).flatten(), np.array(v).flatten(),
                    sigma=selected_data_sigma[k],
                    absolute_sigma=True,
                    method='trf'
                    )
                baseline1[k][i] = fit[0][0]
                baseline1_var[k][i] = fit[1][0, 0] + np.var(v)
            elif v.size > 0:
                baseline1[k][i] = v.values.mean()
                check.loc[i, f'{k}_flag'] = 2
            else:
                baseline1[k][i] = np.nan
                check.loc[i, f'{k}_flag'] = 3
            check.loc[i, f'{k}_cnt'] = v.size
            check.loc[i, f'{k}_win_s'] = date_range0[k]
            check.loc[i, f'{k}_win_e'] = date_range1[k]

    for k, v in baseline1.items():
        v.interpolate('time', inplace=True, limit_direction='both')
    for k, v in baseline1_var.items():
        v.interpolate('time', inplace=True, limit_direction='both')

    baseline2 = {
            k: v.rolling(  # Rolling mean of fitted data
                window2+1, center=True
                ).mean().interpolate(
                    'time', limit_direction='both'
                    )
            for k, v in baseline1.items()
            }
    baseline2_var = {  # Rolling mean of variance
            k: v.rolling(
                window2+1, center=True
                ).mean().interpolate(
                    'time', limit_direction='both'
                    )
            for k, v in baseline1_var.items()
            }
    name_qch4_couple.io.w_nc(
        dim={
            'time': [
                (
                    dates_tHour[dates_tHour.year == year]
                    - pd.to_datetime('1990-01')
                    ).astype('timedelta64[s]'),
                {'size': None},
                {'datatype': np.float64, 'zlib': True},
                {
                    'units': f'seconds since {pd.to_datetime("1990-01")}',
                    'long_name': 'time',
                    'calendar': 'gregorian',
                    }
                ],
            },
        var={
            **{
                k: [
                    np.array(v.loc[f'{year}']).flatten(),
                    {
                        'dimensions': ('time'), 'datatype': np.float64,
                        'zlib': True
                        },
                    {
                        'units': 'ppb' if k == 'chi' else 'per mille',
                        'long_name': (
                            'CH4 mixing ratio' if k == 'chi' else
                            'd13C-CH4' if k == 'd_13C' else
                            'd2H-CH4'
                            )
                        }
                    ]
                for k, v in baseline2.items()
                },
            **{
                f'var_{k}': [
                    np.array(v.loc[f'{year}']).flatten(),
                    {
                        'dimensions': ('time'), 'datatype': np.float64,
                        'zlib': True
                        },
                    {
                        'units': 'ppb' if k == 'chi' else 'per mille',
                        'long_name': (
                            'CH4 mixing ratio_variance' if k == 'chi' else
                            'd13C-CH4_variance' if k == 'd_13C' else
                            'd2H-CH4_variance'
                            )
                        }
                    ]
                for k, v in baseline2_var.items()
                },
            },
        ofile=os.path.join(odir, ofile_fin_nc)
        )
    for k, v in baseline2.items():
        v.to_hdf(os.path.join(odir, ofile_fin), f'{k}')
    for k, v in baseline2_var.items():
        v.to_hdf(os.path.join(odir, ofile_fin), f'var_{k}')

polluted = {
    k:
    {k1: v[obs_filter[k][k1] == 0.] for k1 in filters_fin}
    for k, v in obs[site].items()
    }

'''
intempre = {
        k:
        name_qch4_couple.io.r_intem_base(
            f'baseline/MH_G_ch4{v}_day.txt.gz', 'D'
            )[idx_H[0]:idx_H[-1]]*-1
        for k, v in {'d_13C':'c13', 'd_2H':'h2'}.items()
        } if site == 'mhd' else {}
'''

plot_params = {
    'chi': {
        'ytick': [[1800., 2300.], 100.],
        'lbl': u'$\chi$ (CH$_4$) (ppb)',
        },
    'd_13C': {
        'ytick': [[-48.5, -47.], 0.5],
        'lbl': u'$\delta$ $^{13}$CH$_4$ (\u2030)',
        },
    'd_2H': {
        'ytick': [[-120., -90.], 5.],
        'lbl': u'$\delta$ $^{12}$CH$_3$D (\u2030)',
        }
    }

plt.close('all')
figs = {k: plt.figure(figsize=(6, 3), dpi=300) for k in baseline2}
figs['map'] = plt.figure(figsize=(6, 3), dpi=300)
for fig in figs.values():
    fig.clf()


def plot_baselines(figs, baseline, polluted, xlim, plot_params):
    cmap = matplotlib.cm.get_cmap('plasma_r')
    for k, v in baseline2.items():
        p = polluted[k]
        name_qch4_couple.plot.generic(
            fig=figs[k],
            idata={
                0: [  # filter - domain edge 
                    'line',
                    [p['d_bg2'].index, p['d_bg2'], 'o'],
                    {'c': '#A6DAF3', 'ms': 1.0 if k == 'chi' else 3.0,
                     'mew': 0., 'zorder': 0},
                    ],
                1: [  # filter - south
                    'line',
                    [p['f_south'].index, p['f_south'], 'o'],
                    {'c': cmap(3/11), 'ms': 1.0 if k == 'chi' else 3.0,
                     'mew': 0., 'zorder': 1},
                    ],
                2: [  # filter - eu
                    'line',
                    [p['f_europe'].index, p['f_europe'], 'o'],
                    {'c': cmap(1/11), 'ms': 1.0 if k == 'chi' else 3.0,
                     'mew': 0., 'zorder': 2},
                    ],
                3: [  # filter - local
                    'line',
                    [p['f_local'].index, p['f_local'], 'o'],
                    {'c': cmap(5/11), 'ms': 1.0 if k == 'chi' else 3.0,
                     'mew': 0., 'zorder': 3},
                    ],
                10: [  # first filter
                    'line',
                    [baseline0[k].index, baseline0[k], 'o'],
                    {'c': '#000000', 'ms': 1.0 if k == 'chi' else 3.0,
                     'mew': 0., 'zorder': 10},
                    ],
                11: [  # final
                    'line',
                    [v.index, v, '-'],
                    {'c': '#0000FF', 'lw': 0.5, 'zorder': 11},
                    ],
                12: [  # final lowerbound
                    'line',
                    [v.index, v - baseline2_var[k]**.5, '--'],
                    {'c': '#0000FF', 'lw': 0.5, 'zorder': 12},
                    ],
                13: [  # final upperbound
                    'line',
                    [v.index, v + baseline2_var[k]**.5, '--'],
                    {'c': '#0000FF', 'lw': 0.5, 'zorder': 13},
                    ],
                #20: [  # InTEM
                #    'line',
                #    [intempre[k].index, intempre[k]['Conc'], '-'],
                #    {'c': '#FF0000', 'lw': 0.5, 'zorder': 20},
                #    ] if k in intempre else ['', [], {}],
                },
            xlim=[
                np.datetime64(xlim[0]) - np.timedelta64(1, 'D'),
                np.datetime64(xlim[-1]) + np.timedelta64(1, 'D'),
                ],
            ylim=plot_params[k]['ytick'][0],
            texts=[
                {
                    'x':0.05/6 , 'y':1.55/3 , 's': plot_params[k]['lbl'],
                    'fontsize': 8,
                    'ha': 'left', 'va': 'center', 'rotation': 90
                    },
                ],
            yticks=np.arange(
                plot_params[k]['ytick'][0][0],
                plot_params[k]['ytick'][0][1] + plot_params[k]['ytick'][1],
                plot_params[k]['ytick'][1]
                ),
            loc_plot=[0.6/6.0, 0.2/3.0, 5.1/6.0, 2.75/3.0],
            tick_fontsize=8,
            xtick_params=[True, False, False],
            )


map_filter = np.zeros((nlat, nlon))
map_filter[np.where(filter1['f_south'])] = 3.0
map_filter[np.where(filter1['f_europe'])] = 1.0
map_filter[np.where(filter1['f_local'])] = 5.0
name_qch4_couple.plot.geographical(
    fig=figs['map'],
    idata=map_filter,
    grid=[
        grid_vertex[0],
        grid_vertex[1]
        ],
    texts=[],
    loc_plot=[
        0.05/5.0,
        0.05/3.0,
        4.90/5.0,
        2.90/3.0
        ],
    loc_bar=(),
    bbox=[
        grid_vertex[0][0], grid_vertex[0][-1],
        grid_vertex[1][0], grid_vertex[1][-1]
        ],
    tick_params = {},
    shapefiles=[
        '/home/ec5/name_qch4_couple/shape_files/'
        'CNTR_BN_10M_2016_4326.shp'
        ],
    colour='plasma_r',
    vlim=[1, 12, np.arange(12), np.arange(12)],
    gridlines=[np.arange(-180, 180, 30), np.arange(-90, 90, 30)],
    extend='neither',
    )
plot_baselines(
    figs, baseline2, polluted, [dates_tHour[0], dates_tHour[-1]], plot_params
    )
for k, fig in figs.items():
    fig.savefig(os.path.join(odir, f'baseline_{site}_{k}_{year}.jpeg'))

