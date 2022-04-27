r"""

Baseline Calculation for Hydrogen
@author Filoteea Moldovan
Adapted from Edward Chung

"""
# Standard Library imports
import argparse
import cartopy.crs as ccrs
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
import name_qch4_couple.plot_h2
import name_qch4_couple.region_EU
import name_qch4_couple.routines
import name_qch4_couple.util

# Local imports
import routines
import chem_co

# =============================================================================
#   Settings
# =============================================================================
# Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument("-site", required=True)
parser.add_argument("-species", required=True)
parser.add_argument("-year", required=True, type=int)
parser.add_argument("-window1", required=True, type=int)
parser.add_argument("-window2", required=True, type=int)
parser.add_argument("-force", required=True, type=int)
parser.add_argument("-odir", required=True)
args = parser.parse_args()

site = args.site
year = args.year
species = args.species
window1 = args.window1
window2 = args.window2
force_compute = bool(args.force)
odir = args.odir

#mode = 'all'
site_ref = 'mhd'  # HDF reference key
p = 101325  # Pa
T = 288.15  # K

ofile_con = f'condition-{site}-{species}-{year}.nc'
ofile_fil = f'filtered-{site}-{species}-{year}.nc'
ofile_fin = f'baseline-{site}-{species}-{year}.nc'

long_names = OrderedDict()
locations = OrderedDict()
with open(f'inputs/baseline/{site}.json', 'r') as f:
    st_info = json.load(f)
    long_names[site] = st_info['long_name']
    locations[site] = st_info['location']

site1 = site.replace('_', '-')

date_nodash = 'REPLACE'
if species == 'h2':
    #from chem_co import read_Q, read_obs, read_baseline
    var_name="chi_H2"
    Dfile = (
        f'inputs/baseline/footprints_mhd/'
        f'{site1}_UKV_EUROPE_{date_nodash}.nc'
        )
    Q2obs = 1.e9 / 2.016  # M_H2 (g mol-1) - IUPAC
    ylabel = u'$\chi$ H$_{2}$ (nmol mol$^{-1}$)'
    ylim = [460., 600.]
    yticks = np.arange(460., 600., 30.)
    var_long_name = 'mole_fraction_of_hydrogen'
    var_units = 'nmol mol-1'

# =============================================================================
#   Pre-processing
# =============================================================================
print(f'Initialising')
print(f'    site = {site}')
print(f'    year = {year}')

# Dates
dt1 = pd.to_timedelta(window1//2, 'H')
dt2 = pd.to_timedelta(window2//2, 'H')
dt_large = max([dt1, dt2]) * 3

dates_tHour, dates_tDay = (
    pd.date_range(
        pd.to_datetime(f'{year}') - dt_large,
        pd.to_datetime(f'{int(year) + 1}') + dt_large,
        freq='1H',
        closed='left'
        )
    for freq in ['1H', '1D']
    )
dates_tMonth =  pd.date_range(
    (pd.to_datetime(f'{year}') - dt_large).strftime('%Y-%m'),
    (pd.to_datetime(f'{int(year) + 1}') + dt_large).strftime('%Y-%m'),
    freq='1MS',
    #closed='left'
    )

# Grid
grid_info = routines.define_grid()
dlat = grid_info['dlat']
nlat = grid_info['nlat']
nlon = grid_info['nlon']
area = grid_info['area']
grid_centre = grid_info['grid_centre']
grid_vertex = grid_info['grid_vertex']

nlat_odd = bool(nlat % 2)
nlon_odd = bool(nlon % 2)
nlat_half = nlat // 2
nlon_half = nlon // 2

# Fooprints
#print('Read Footprints')
def read_D(Dfile):
    h1 = 3.e3
    h2 = 8.e3
    with xr.open_dataset(Dfile) as ds_read:
        with ds_read.load() as Din:
            # time
            nt = Din.time.size
            # latlon
            nlat = Din.lat.size
            nlat_odd = bool(nlat % 2)
            nlat_half = nlat // 2
            nlon = Din.lon.size
            nlon_odd = bool(nlon % 2)
            nlon_half = nlon // 2
            # height
            lt = Din.height.values < h1
            ut = (h1 <= Din.height.values) & (Din.height.values < h2)
            st = h2 <= Din.height.values
            # Footprint (lat, lon, time)
            D = Din.fp.transpose('time', 'lat', 'lon').values
            # End locations (height, lat/lon, time)
            end = np.zeros((nt, 17))
            endn = Din.particle_locations_n.transpose('time', 'height', 'lon'
                    ).values
            ends = Din.particle_locations_s.transpose('time', 'height', 'lon'
                    ).values
            ende = Din.particle_locations_e.transpose('time', 'height', 'lat'
                    ).values
            endw = Din.particle_locations_w.transpose('time', 'height', 'lat'
                    ).values
            end[:, 0] += endn[:, lt, -nlon_half:].sum((1, 2))
            end[:, 1] += ende[:, lt, :+nlat_half].sum((1, 2))
            end[:, 2] += ende[:, lt, -nlat_half:].sum((1, 2))
            end[:, 3] += ends[:, lt, :+nlon_half].sum((1, 2))
            end[:, 4] += ends[:, lt, -nlon_half:].sum((1, 2))
            end[:, 5] += endw[:, lt, :+nlat_half].sum((1, 2))
            end[:, 6] += endw[:, lt, -nlat_half:].sum((1, 2))
            end[:, 7] += endn[:, lt, :+nlon_half].sum((1, 2))
            end[:, 8] += endn[:, ut, -nlon_half:].sum((1, 2))
            end[:, 9] += ende[:, ut, :+nlat_half].sum((1, 2))
            end[:, 10] += ende[:, ut, -nlat_half:].sum((1, 2))
            end[:, 11] += ends[:, ut, :+nlon_half].sum((1, 2))
            end[:, 12] += ends[:, ut, -nlon_half:].sum((1, 2))
            end[:, 13] += endw[:, ut, :+nlat_half].sum((1, 2))
            end[:, 14] += endw[:, ut, -nlat_half:].sum((1, 2))
            end[:, 15] += endn[:, ut, :+nlon_half].sum((1, 2))
            end[:, 16] += (
                endn[:, st].sum((1, 2))
                + ende[:, st].sum((1, 2))
                + ends[:, st].sum((1, 2))
                + endw[:, st].sum((1, 2))
                )
            if nlon_odd:
                end[:, 0] += endn[:, lt, nlon_half].sum((1,)) / 2
                end[:, 1] += ende[:, lt, nlat_half].sum((1,)) / 2
                end[:, 2] += ende[:, lt, nlat_half].sum((1,)) / 2
                end[:, 3] += ends[:, lt, nlon_half].sum((1,)) / 2
                end[:, 4] += ends[:, lt, nlon_half].sum((1,)) / 2
                end[:, 5] += endw[:, lt, nlat_half].sum((1,)) / 2
                end[:, 6] += endw[:, lt, nlat_half].sum((1,)) / 2
                end[:, 7] += endn[:, lt, nlon_half].sum((1,)) / 2
                end[:, 8] += endn[:, ut, nlon_half].sum((1,)) / 2
                end[:, 9] += ende[:, ut, nlat_half].sum((1,)) / 2
                end[:, 10] += ende[:, ut, nlat_half].sum((1,)) / 2
                end[:, 11] += ends[:, ut, nlon_half].sum((1,)) / 2
                end[:, 12] += ends[:, ut, nlon_half].sum((1,)) / 2
                end[:, 13] += endw[:, ut, nlat_half].sum((1,)) / 2
                end[:, 14] += endw[:, ut, nlat_half].sum((1,)) / 2
                end[:, 15] += endn[:, ut, nlon_half].sum((1,)) / 2
    return D, end

# Observations

def r_decc(fpath):

    odata = pd.read_csv(
                fpath,
                usecols=lambda x: x.lower() in ['time', 'h2_ppb'],
                index_col=['time'],
                skipinitialspace=True,
                parse_dates=['time']
                ).dropna()
#    odata = odata.dropna()
    odata.columns = odata.columns.str.lower()
    return odata

def read_obs(timestamps, site, resample=False):
    date = timestamps[0].strftime('%Y-%m')
    t0 = timestamps[0].strftime('%Y-%m-%d %H')
    t1 = timestamps[-1].strftime('%Y-%m-%d %H')
    if site == 'WAO':
        ifile = 'inputs/obs/WAO_H2_oct2021.csv'
        col_or_no = 'h2_ppb'
        sigma_col_or_no = 0.2
    elif site == 'MHD_10magl':
        ifile = 'inputs/baseline/MHD_2018.csv'
        col_or_no = 'h2_ppb'
        sigma_col_or_no = 0.2
    else:
        ifile = False
        col_or_no = np.nan
        sigma_col_or_no = np.nan
    if ifile:
        all_obs_raw = r_decc(ifile).sort_index().loc[t0:t1]
        obs_raw = all_obs_raw[col_or_no]
        sigma_obs_raw = (all_obs_raw[sigma_col_or_no]
                         if isinstance(sigma_col_or_no, str) else
                         pd.Series(sigma_col_or_no, index=all_obs_raw.index))
    if isinstance(col_or_no, str):
        obs = (obs_raw
               if resample is False else
               obs_raw.resample('1H').mean().reindex(timestamps))
    else:
        obs = pd.Series(col_or_no, index=timestamps)
    if isinstance(sigma_col_or_no, str) or isinstance(col_or_no, str):
        sigma_obs = (
            sigma_obs_raw
            if resample is False else
            sigma_obs_raw.resample('1H').apply(
                lambda x: np.sum(x**2)).reindex(timestamps))
    else:
        sigma_obs = pd.Series(sigma_col_or_no, index=timestamps)
    return obs, sigma_obs 


obs, sigma_obs = read_obs(dates_tHour, site)

# Filters - Footprint (should not exceed)
filter1 = {
    k: np.zeros((nlat, nlon)) for k in ['f_local', 'f_europe', 'f_south']
    }

## Local
radius = 2
st_idx = [np.abs(grid_centre[i] - l).argmin() for i, l in enumerate(locations[site])]
filter1['f_local'][
    st_idx[1]-radius:st_idx[1]+radius+1,
    st_idx[0]-radius:st_idx[0]+radius+1
    ] = 1.

## Europe
filter1['f_europe'][
    np.load('europe.npy', 'r')
    ] = 1.

## Southern
lat_south = 30.
filter1['f_south'][grid_centre[1] + 0.5*dlat <= lat_south] = 1.

## Population
# filter1['f_pop'] = 

## Thresholds
filter_threshold = {k: 1. for k in filter1}
filter_threshold['f_local'] = 10
filter_threshold['f_europe'] = 5
filter_threshold['f_south'] = 2

# Dilution sensitivity limit
#dsl = 3.4 * 1.e6 / p * 8.314 * T
dsl = 3.4 / p * 8.314 * T

# Filters - Domain borders (should exceed)
filter2 = {
    'p_nw_all': [0, 5, 6, 7, 8, 13, 14, 15],  # default
    }
## Thresholds
filter_threshold['p_nw_all'] = 0.8

# =============================================================================
#   Main Process
# =============================================================================
print("Processing Main Routine")
filters = list(filter1) + list(filter2)
flag_con = (ofile_con
            and os.path.exists(os.path.join(odir, ofile_con))
            and not force_compute)
if flag_con:
    try:
        with xr.open_dataset(os.path.join(odir, ofile_con)) as ds_read:
            with ds_read.load() as ds:
                if any([i not in ds.variables for i in filter_threshold]):
                    flag_con = False
    except:
        flag_con = False
                    
# Get conditions
if flag_con:
    print("    Using pre-existing filters")
    with xr.open_dataset(os.path.join(odir, ofile_con)):
        with ds_read.load() as ds:
            condition = pd.concat([
                ds[i].to_series() for i in filter_threshold
                ], axis=1)
            idx_H = condition.index
else:
    print("    Applying filters")
    idx_H = dates_tHour 
    condition = pd.DataFrame(0, index=idx_H, columns=filters)
    for month in dates_tMonth:
        print(f"        {month}")
        month_idx = month.strftime('%Y-%m')
        try:
            hoursinmonth = pd.date_range(
                month,
                month + pd.offsets.MonthBegin(1),
                freq='1H',
                closed='left'
                )
            Dfile_in = re.sub(date_nodash, month.strftime('%Y%m'), Dfile)
            D, end_loc = (i[hoursinmonth.isin(dates_tHour)] for i in read_D(Dfile_in))
            condition.loc[month_idx, filter1] = pd.DataFrame(
                {k: (D * v).sum((1, 2)) for k, v in filter1.items()},
                index=condition.loc[month_idx].index
                )
            condition.loc[month_idx, filter2] = pd.DataFrame(
                {k: end_loc[:, v].sum((1,)) for k, v in filter2.items()},
                index=condition.loc[month_idx].index
                )
        except:
            condition.loc[month, filter1] = np.inf
            condition.loc[month, filter2] = -np.inf
    name_qch4_couple.io.w_nc(
        dim={
            'time': [
                (
                    dates_tHour
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
                    condition[k].values,
                    {
                        'dimensions': ('time'), 'datatype': np.float64,
                        'zlib': True
                        },
                    {
                        'units': var_units,
                        'long_name': var_long_name
                        }
                    ]
                for k in filter_threshold
                },
            },
        ofile=os.path.join(odir, ofile_con)
        )

# Check if baseline condition
filtered = pd.DataFrame(0, index=idx_H, columns=filters)
for k in filter1:
    threshold = filter_threshold[k]
    filtered.loc[condition[k] / dsl <= threshold, k] = 1.
for k in filter2:
    threshold = filter_threshold[k]
    filtered.loc[condition[k] >= threshold, k] = 1.

# Baseline
filters_fin = ['f_local', 'f_europe', 'f_south', 'p_nw_all']
oidx = obs.index
nidx = oidx.floor('H')
obs_filter = filtered.reindex(nidx).reset_index(drop=True)
obs_filter.index = oidx
baseline0 = obs[obs_filter[filters_fin].min(1) == 1.]
baseline0_sigma = sigma_obs[obs_filter[filters_fin].min(1) == 1.]

#baseline0 = baseline0.loc[baseline0.notnull()]
#baseline0_sigma = baseline0_sigma.loc[baseline0.notnull()]


if not force_compute and os.path.exists(os.path.join(odir, ofile_fin)):
    print("    Reading baselines")
    with xr.open_dataset(os.path.join(odir, ofile_fin)) as ds_read:
        with ds_read.load() as ds:
            baseline2 = ds[var_name].to_series()
            baseline2_var = ds[f'var_{var_name}'].to_series()
else:
    print("    Processing baselines")
    warnings.simplefilter('ignore', np.RankWarning)
    check = pd.DataFrame(0., index=dates_tHour, columns=[])

    meas_on_time = dates_tHour.intersection(baseline0.index)

    window1bef = dates_tHour.to_series().apply(
        lambda x: baseline0.loc[x - dt1//2 : x].size
        )
    window1aft = dates_tHour.to_series().apply(
        lambda x: baseline0.loc[x : x + dt1//2].size
        )
    window1minhalfwidth = dates_tHour.to_series().apply(
        lambda x: 2 if x in meas_on_time else 1
        )
    window1bef_more = window1minhalfwidth - window1bef
    window1aft_more = window1minhalfwidth - window1aft

    window1start = dates_tHour.to_series().apply(lambda x: (
            x - dt1
            if window1bef.loc[x] > window1minhalfwidth.loc[x] else
            baseline0.loc[:x - dt1].iloc[-window1bef_more[x]:].first_valid_index()
            if baseline0.loc[:x - dt1].iloc[-window1bef_more[x]:].size else
            dates_tHour[0]
            ))
    window1end = dates_tHour.to_series().apply(
        lambda x: (
            x + dt1
            if window1aft.loc[x] > window1minhalfwidth.loc[x] else
            baseline0.loc[x + dt1:].iloc[:window1aft_more[x]].last_valid_index()
            if baseline0.loc[x + dt1:].iloc[:window1aft_more[x]].size else
            dates_tHour[-1]
            ))

    def fit_func(value, sigma, time):
        nfront = value.loc[:time].size
        nback = value.loc[time:].size
        if value.size >= 2 and all([nfront, nback]):
            order = min(value.size, 2)
            t2H = (value.index - time).to_series() / pd.Timedelta(hours=1)
            fit = scipy.optimize.curve_fit(
                (
                    (lambda x, a, b, c: a + b*x + c*x**2)
                    if order == 2 else
                    (lambda x, a, b: a + b*x)
                    ),
                t2H.values, value.values,
                sigma=sigma.values,
                absolute_sigma=True,
                method='trf'
                )
            intercept = fit[0][0]
            variance = fit[1][0, 0] + np.var(value)
        elif value.size:
            intercept = value.values[0]
            variance = sigma.values[0]
        else:
            intercept = np.nan
            variance = np.nan
        return pd.Series([intercept, variance], index=['intercept', 'variance'])

    baseline1 = dates_tHour.to_series().apply(lambda x:
            fit_func(baseline0.loc[window1start[x]:window1end[x]],
                     baseline0_sigma.loc[window1start[x]:window1end[x]],
                     x))

    baseline1.interpolate('time', inplace=True, limit_direction='both')

    baseline2 = baseline1['intercept'].rolling(
        window2+1, center=True
        ).mean().interpolate('time', limit_direction='both')
    baseline2_var = baseline1['variance'].rolling(
        window2+1, center=True
        ).mean().interpolate('time', limit_direction='both')
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
            var_name: [
                np.array(baseline2.loc[f'{year}']).flatten(),
                {
                    'dimensions': ('time'), 'datatype': np.float64,
                    'zlib': True
                    },
                {
                    'units': var_units,
                    'long_name': var_long_name
                    }
                ],
            f'var_{var_name}': [
                np.array(baseline2_var.loc[f'{year}']).flatten(),
                {
                    'dimensions': ('time'), 'datatype': np.float64,
                    'zlib': True
                    },
                {
                    'units': var_units,
                    'long_name': f'{var_long_name}_variance'
                    }
                ],
            },
        ofile=os.path.join(odir, ofile_fin)
        )

polluted = {k: obs[obs_filter[k] == 0.] for k in filters_fin}

# Plot
figs = {}
axs = {}
pobjs = {}

cmap = matplotlib.cm.get_cmap('plasma_r')

colours = {
    'f_local': cmap(5/11),
    'f_europe': cmap(1/11),
    'f_south': cmap(3/11),
    'p_nw_all': '#A6DAF3',
    'background': '#000000',
    'final': '#0000FF'
    }
zorder = {
    'f_local': 4,
    'f_europe': 3,
    'f_south': 2,
    'p_nw_all': 1,
    'background': 5,
    'final': 6
    }
fig_param = {
    'mw': 6, 'mh': 5,
    'mpw': 5.3, 'mph': 4.5,
    'mgap': 0.05,
    'mlmargin': 0.7, 'mbmargin': 0.5,
    'ylblx': 0.05, 'ylbly': 1.5,  # left, centre aligned
    'fontsize': 8,
    }
plt.close('all')

figs['main'] = plt.figure(figsize=(fig_param['mw'], fig_param['mh']), dpi=300)
axs['main'] = {}
pobjs['main'] = {}
figs['main'].clf()

name_qch4_couple.plot_h2.generic3(
    fig=figs['main'],
    axs=axs['main'],
    pobjs=pobjs['main'],
    new_axs={
        'date1': [
            dict(
                rect=[
                    (1 * fig_param['mlmargin']
                        + fig_param['mgap'])
                    / fig_param['mw'],
                    (fig_param['mh']
                        - fig_param['mph']
                        + fig_param['mgap'])
                    / fig_param['mh'],
                    (fig_param['mpw']
                        - 2*fig_param['mgap'])
                    / fig_param['mw'],
                    (fig_param['mph']
                        - 2*fig_param['mgap'])
                    / fig_param['mh']
                    ],
                label='date1',
                projection=None
                ),
            {
                "set_yticks": [[yticks], {}],
                "set_xlim": [[
                    dates_tHour[dates_tHour.year == year][0] - pd.to_timedelta('1 day'),
                    dates_tHour[dates_tHour.year == year][-1] + pd.to_timedelta('1 day')
                    ], {}],
                "set_ylim": [[ylim], {}],
                "tick_params": [[], dict(
                    axis='both', which='major', direction='in',
                    labelsize=fig_param['fontsize'],
                    left=True, bottom=True,
                    right=False, top=False,
                    labelleft=True, labelbottom=True,
                    labelright=False, labeltop=False,
                    )],
                "xaxis.set_major_locator": [
                    [mdates.DayLocator(bymonthday=1)], {}
                    ],
                "xaxis.set_major_formatter": [
                    [mdates.DateFormatter('%Y-%m-%d')], {}
                    ],
                },
            dict(
                patch_alpha=0.0
                )
            ]
        },
    new_pobjs={
        **{
            k: [
                'date1', 'plot',
                [polluted[k].index, np.array(polluted[k]), 'o'],
                {'c': colours[k], 'ms': 1., 'mew': 0., 'zorder': zorder[k],
                 'label': k}
                ]
            for k in polluted
            },
        'background': [
            'date1', 'plot', [baseline0.index, baseline0.values, 'o'],
            {'c': colours['background'], 'ms': 1., 'mew': 0.,
             'zorder': zorder['background'],
             'label': 'background'}
            ],
        'final': [
            'date1', 'plot', [baseline2.index, baseline2.values, '-'],
            {'c': colours['final'], 'ms': 1., 'mew': 0.,
             'zorder': zorder['final'],
             'label': 'final'}
            ],
        'final_l': [
            'date1', 'plot',
            [baseline2.index, (baseline2 - baseline2_var**.5).values, '--'],
            {'c': colours['final'], 'ms': 1., 'mew': 0.,
             'zorder': zorder['final'],
             'label': 'final (lower)'}
            ],
        'final_u': [
            'date1', 'plot',
            [baseline2.index, (baseline2 + baseline2_var**.5).values, '--'],
            {'c': colours['final'], 'ms': 1., 'mew': 0.,
             'zorder': zorder['final'],
             'label': 'final (upper)'}
            ],
        'legend':[
            'date1', 'legend', [],
            dict(
                loc='upper right',
                numpoints=2, fontsize=fig_param['fontsize'], ncol=3,
                markerscale=5.0/3.5, handletextpad=0.2, columnspacing=1.0,
                borderpad=0.2, borderaxespad=0.2
                )
            ]
        },
    texts=[
        {
            'x': fig_param['mgap'] / fig_param['mw'],
            'y': (fig_param['mh']
                    - 1/2*fig_param['mph'])
                / fig_param['mh'],
            's': ylabel,
            'ha': 'left', 'va': 'center',
            'size': fig_param['fontsize'], 'rotation': 90
            }
        ],
    legend_params=[
        [],
        [],
        {}
        ]
    )
for l in axs['main']['date1'].get_xticklabels():
    l.set_ha("right")
    l.set_rotation(30)
figs['main'].savefig(os.path.join(odir, f'baseline-{site}-{species}-{year}mhd_magl.png'))

map_filter = np.zeros((nlat, nlon))
map_filter[np.where(filter1['f_south'])] = 3.0
map_filter[np.where(filter1['f_europe'])] = 1.0
map_filter[np.where(filter1['f_local'])] = 5.0
fig_param['mh'] = 3
figs['map'] = plt.figure(figsize=(fig_param['mw'], fig_param['mh']), dpi=300)
figs['map'].clf()

'''
name_qch4_couple.plot_h2.geographical(
    fig=figs['map'],
    idata=map_filter,
    grid=[
        grid_vertex[0],
        grid_vertex[1]
        ],
    texts=[],
    loc_plot=[
        fig_param['mgap'] / fig_param['mw'],
        fig_param['mgap'] / fig_param['mh'],
        (fig_param['mw'] - 2*fig_param['mgap']) / fig_param['mw'],
        (fig_param['mh'] - 2*fig_param['mgap']) / fig_param['mh'],
        ],
    loc_bar=(),
    bbox=[
        grid_vertex[0][0], grid_vertex[0][-1],
        grid_vertex[1][0], grid_vertex[1][-1]
        ],
    tick_params = {},
    shapefiles=['C:/Users/filot/Desktop/YEAR_4/Dissertation/Ed_new_script/scripts/shape_files/CNTR_BN_10M_2016_4326.shp'],
    colour='plasma_r',
    vlim=[1, 12, np.arange(12), np.arange(12)],
    gridlines=[np.arange(-180, 180, 30), np.arange(-90, 90, 30)],
    extend='neither',
    )

name_qch4_couple.plot_h2.geographical(
    fig=figs['map'],
    idata=map_filter,
    grid=[
        grid_vertex[0],
        grid_vertex[1]
        ],
    texts=[],
    loc_plot=[
        0.3 / fig_param['mw'],
        0.9 / fig_param['mh'],
        3.5 / fig_param['mw'],
        2.0 / fig_param['mh'],
        ],
    loc_bar=(),
    bbox=[
        -12, 6, 49, 61 
        ],
    tick_params = {},
    shapefiles=['C:/Users/filot/Desktop/YEAR_4/Dissertation/Ed_new_script/scripts/shape_files/CNTR_BN_10M_2016_4326.shp'],
    colour='plasma_r',
    vlim=[1, 12, np.arange(12), np.arange(12)],
    gridlines=[np.arange(-180, 180, 30), np.arange(-90, 90, 30)],
    extend='neither',
    )
figs['map'].savefig(os.path.join(odir, f'baseline-{site}-{species}-{year}-map.png'))
'''
