#!/p/system/packages/anaconda/5.0.0_py3/bin/python3.6
# run with Python 3.6

import argparse
import os
import fnmatch
import pandas as pd
import xarray as xr
from pathlib import Path

parser = argparse.ArgumentParser(description='Extract climate data for local lakes sector.')

parser.add_argument('-f', '--lakes-file', dest='lakes_file',
                    default='lakes.csv',
                    help='CSV file with lake long name, specifier, latitude and longitude.')
parser.add_argument('-p', '--phase', dest='phase', required=True,
                    default='ISIMIP3a',
                    help='ISIMIP phase.')
parser.add_argument('-t', '--tier', dest='tier', required=True,
                    default='InputData',
                    help='ISIMIP Input data tier, e.g. [InputData|SecondaryInputData]')
parser.add_argument('-d', '--datatype', dest='datatype', required=True,
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

path = args.basedir + '/' + args.phase + '/' + args.tier + '/climate/atmosphere/' \
       + args.datatype + '/global/daily/' + '/' + args.climforcing + '/' + args.model

outpath = Path(args.outdir)

lakes_df = pd.read_csv(args.lakes_file)


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
    tas_ds = xr.load_dataset(path + '/' + tas_file)
    hurs_ds = xr.load_dataset(path + '/' + hurs_file)
    pr_ds = xr.load_dataset(path + '/' + pr_file)
    rsds_ds = xr.load_dataset(path + '/' + rsds_file)
    rlds_ds = xr.load_dataset(path + '/' + rlds_file)
    ps_ds = xr.load_dataset(path + '/' + ps_file)
    sfcwind_ds = xr.load_dataset(path + '/' + sfcwind_file)

    # iterate over lakes
    for index, lake in lakes_df.iterrows():
        print('  {1}, lat:{2} lon:{3}'.format(*lake))

        outfile = args.model.lower() + '_' + args.climforcing + '_' + args.datatype + '_' + \
            lake[1].replace(' ', '-').lower() + '_daily_' + first_year + '_' + last_year + '.csv'

        # read actual data from NetCDF
        print('   read tas ...')
        tas_lake = tas_ds.sel(lat=lake[2], lon=lake[3], method='nearest', drop=True)
        print('   read hurs ...')
        hurs_lake = hurs_ds.sel(lat=lake[2], lon=lake[3], method='nearest', drop=True)
        print('   read pr ...')
        pr_lake = pr_ds.sel(lat=lake[2], lon=lake[3], method='nearest', drop=True)
        print('   read rsds ...')
        rsds_lake = rsds_ds.sel(lat=lake[2], lon=lake[3], method='nearest', drop=True)
        print('   read rlds ...')
        rlds_lake = rlds_ds.sel(lat=lake[2], lon=lake[3], method='nearest', drop=True)
        print('   read ps ...')
        ps_lake = ps_ds.sel(lat=lake[2], lon=lake[3], method='nearest', drop=True)
        print('   read sfcwind ...')
        sfcwind_lake = sfcwind_ds.sel(lat=lake[2], lon=lake[3], method='nearest', drop=True)

        # write csv files per lake
        print('   write data ...')

        lake_data = tas_lake.to_dataframe()
        lake_data['hurs'] = hurs_lake['hurs']
        lake_data['pr'] = pr_lake['pr']
        lake_data['rsds'] = rsds_lake['rsds']
        lake_data['rlds'] = rlds_lake['rlds']
        lake_data['ps'] = ps_lake['ps']
        lake_data['sfcwind'] = sfcwind_lake['sfcwind']

        if period.split('_')[0] == first_year:
            lake_data.to_csv(outpath / outfile)
        else:
            lake_data.to_csv(outpath / outfile, mode='a', header=False)

    tas_ds.close()
    hurs_ds.close()
    pr_ds.close()
    rsds_ds.close()
    rlds_ds.close()
    ps_ds.close()
    sfcwind_ds.close()
