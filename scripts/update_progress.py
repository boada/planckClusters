from glob import glob
from numpy import where, inf
import pandas as pd

def update_sheet(spreadsheet, file_list, mode='Mosaics'):
    ''' this updates the spreadsheet with the status of whatever the "mode" is.
    So if mode=mosaics then it updates the status of the mosaics. If it is
    mode=catalogs then it updates the status of the catalogs. In either case
    you have to pass it a file list of the appropriate files. It will fail
    otherwise.

    '''

    # remove any exisiting data
    spreadsheet[mode] = ''

    # set the sort order
    order = ['g', 'r', 'i', 'z', 'K']

    for f in file_list:
        file = f.split('/')[-1]

        if mode == 'Mosaics':
            # remove extention
            file = file[:-5]
            if 'weight' in file:
                continue
            if 'SEG' in file:
                continue
            if 'BACK' in file:
                continue
            band = file[-1]
            obj = file[:-1]

            if 'PSZ' in obj:
                try:
                    indx = where(spreadsheet['Name'] == obj)[0][0]
                    spreadsheet[mode][indx] = spreadsheet[mode][indx] + band

                except IndexError:
                    pass

        elif mode == 'Catalogs':
            # remove extension
            file = file[:-4]
            if 'cal' in file:
                # remove the calibrated portion
                file = file[:-4]
                band = file[-1]
                obj = file[:-1]
            elif 'bpz' in file:
                file = file[:-4]
                band = file[-1]
                obj = file[:-1]
            else:
                band = file[-1]
                obj = file[:-1]

            print(obj)
            if 'PSZ' in obj:
                try:
                    indx = where(spreadsheet['Name'] == obj)[0][0]
                    spreadsheet[mode][indx] = spreadsheet[mode][indx] + band

                except IndexError:
                    pass

        elif mode == 'Calibrated':
            # remove extension
            file = file[:-8]
            # remove the calibrated portion
            band = file[-1]
            obj = file[:-1]
            if 'PSZ' in obj:
                try:
                    indx = where(spreadsheet['Name'] == obj)[0][0]
                    spreadsheet[mode][indx] = spreadsheet[mode][indx] + band

                except IndexError:
                    pass

        elif mode == 'Photo-z':
            # remove extension
            file = file[:-4]
            obj = file
            print(obj)
            if 'PSZ' in obj:
                try:
                    indx = where(spreadsheet['Name'] == obj)[0][0]
                    spreadsheet[mode][indx] = 'True'
                except IndexError:
                    pass

        else:
            print('Mode not understood')
            return

    if mode not in ['Photo-z']:
        for indx in range(len(spreadsheet)):
            spreadsheet[mode][indx] = ''.join(sorted(spreadsheet[mode][indx],
                                         key=lambda x: order.index(x) if x in
                                         order else inf))

    return spreadsheet


spreadsheet = pd.read_csv('/home/boada/Projects/planckClusters/'
                'catalogs/PSZ2_unconfirmed_catalog - Master.csv')

mosaics = glob('./proc2/**/*.fits', recursive=True)
catalogs = glob('./proc2/**/*.cat', recursive=True)
calibrated = glob('./proc2/**/*_cal.cat', recursive=True)
photz = glob('./proc2/**/*.bpz', recursive=True)

spreadsheet = update_sheet(spreadsheet, mosaics, mode='Mosaics')
spreadsheet = update_sheet(spreadsheet, catalogs, mode='Catalogs')
spreadsheet = update_sheet(spreadsheet, calibrated, mode='Calibrated')
spreadsheet = update_sheet(spreadsheet, photz, mode='Photo-z')

spreadsheet.to_csv('updated.csv')

