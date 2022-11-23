#!/p/system/packages/anaconda/5.0.0_py3/bin/python3.6
# run with Python 3.6

from netCDF4 import Dataset, num2date
import argparse
import os
import fnmatch
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser(description='Extract climate data for local lakes sector.')

parser.add_argument('-p', '--phase', dest='phase', required=True,
                    default='ISIMIP3a',
                    help='ISIMIP phase.')
parser.add_argument('-t', '--datatype', dest='datatype', required=True,
                    help='Input data types, e.g. 3a: [obsclim|counterclim], 3b: [bias-corrected]')
parser.add_argument('-m', '--model', dest='model', required=True,
                    help='Input model to process')
parser.add_argument('-c', '--climate-forcing', dest='climforcing', required=True,
                    help='climate forcing model to process, e.g [picontrol|historical|ssp126|ssp370|ssp585].')
parser.add_argument('-b', '--base', dest='basedir',
                    default='/p/projects/isimip/isimip',
                    help='base directory')
parser.add_argument('-o', '--out', dest='outdir',
                    default='/p/tmp/buechner/extract_lakes_python',
                    help='base directory to post-processed data')


args = parser.parse_args()

path = args.basedir + '/' + args.phase + '/InputData/climate/atmosphere/' \
       + args.datatype + '/global/daily/' + '/' + args.climforcing + '/' + args.model

outpath = Path(args.outdir)

lakes_df = pd.read_csv('lakes.csv')


# return available time periods of our daily data chunked to decades
def get_periods(path):
    files = sorted(fnmatch.filter(os.listdir(path), '*_tas_*'))
    periods = [file.split('daily_')[1].split('.nc')[0] for file in files]
    return periods


first_year = get_periods(path)[0].split(sep='_')[0]
last_year = get_periods(path)[-1].split(sep='_')[1]

for period in get_periods(path):
    print(' Period: ', period)

    # get file names for current period
    tas_file = [file for file in fnmatch.filter(os.listdir(path), '*_tas_*') if period in file][0]
    hurs_file = [file for file in fnmatch.filter(os.listdir(path), '*_hurs_*') if period in file][0]
    pr_file = [file for file in fnmatch.filter(os.listdir(path), '*_pr_*') if period in file][0]
    rsds_file = [file for file in fnmatch.filter(os.listdir(path), '*_rsds_*') if period in file][0]
    rlds_file = [file for file in fnmatch.filter(os.listdir(path), '*_rlds_*') if period in file][0]
    ps_file = [file for file in fnmatch.filter(os.listdir(path), '*_ps_*') if period in file][0]
    sfcwind_file = [file for file in fnmatch.filter(os.listdir(path), '*_sfcwind_*') if period in file][0]

    # open datasets
    tas_ds = Dataset(path + '/' + tas_file)
    hurs_ds = Dataset(path + '/' + hurs_file)
    pr_ds = Dataset(path + '/' + pr_file)
    rsds_ds = Dataset(path + '/' + rsds_file)
    rlds_ds = Dataset(path + '/' + rlds_file)
    ps_ds = Dataset(path + '/' + ps_file)
    sfcwind_ds = Dataset(path + '/' + sfcwind_file)

    # get dates
    times = tas_ds.variables['time']
    time_units = times.units
    time_calendar = times.calendar
    dates = num2date(times[:], units=time_units, calendar=time_calendar)

    # iterate over lakes
    for index, lake in lakes_df.iterrows():
        print('{0} (lat:{1} lon:{2})'.format(*lake))

        outfile = args.model.lower() + '_' + args.climforcing + '_' + args.datatype + '_' + \
            lake[0].replace(' ', '-').lower() + '_daily_' + first_year + '_' + last_year + '.csv'

        idx_lat = round(2 * (89.75 - float(lake[1])))
        idx_lon = round(2 * (float(lake[2]) + 179.75))

        # read actual data from NetCDF
        print('   read tas ...')
        tas_lake = tas_ds['tas'][:, idx_lat, idx_lon]
        print('   read hurs ...')
        hurs_lake = hurs_ds['hurs'][:, idx_lat, idx_lon]
        print('   read pr ...')
        pr_lake = pr_ds['pr'][:, idx_lat, idx_lon]
        print('   read rsds ...')
        rsds_lake = rsds_ds['rsds'][:, idx_lat, idx_lon]
        print('   read rlds ...')
        rlds_lake = rlds_ds['rlds'][:, idx_lat, idx_lon]
        print('   read ps ...')
        ps_lake = ps_ds['ps'][:, idx_lat, idx_lon]
        print('   read sfcwind ...')
        sfcwind_lake = sfcwind_ds['sfcwind'][:, idx_lat, idx_lon]

        lake_data = pd.DataFrame({
            'date': dates[:],
            'tas': tas_lake[:],
            'hurs': hurs_lake[:],
            'pr': pr_lake[:],
            'rsds': rsds_lake[:],
            'rlds': rlds_lake[:],
            'ps': ps_lake[:],
            'sfcwind': sfcwind_lake[:]
        })

        # write csv files per lake
        print('   write data ...')
        if period.split('_')[0] == first_year:
            lake_data.to_csv(outpath / outfile, encoding='utf-8', float_format='%.1f', index=False)
        else:
            lake_data.to_csv(outpath / outfile, encoding='utf-8', float_format='%.1f', index=False, mode='a', header=False)

    tas_ds.close()
    hurs_ds.close()
    pr_ds.close()
    rsds_ds.close()
    rlds_ds.close()
    ps_ds.close()
    sfcwind_ds.close()
