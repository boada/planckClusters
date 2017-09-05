import shutil
from glob import glob
import os
import calendar
from astropy.io import fits

''' This file takes a bunch of downloaded VO data and sorts it by the type of
of observation and date. It keeps the resampled/sacked/etc files together. It
should only be run on a big directory of unsorted VO data.

'''

files = glob('*.fz')
for f in files:
    with fits.open(f) as oimg:
        try:
            date = oimg[0].header['Date-obs']
            inst = oimg[0].header['INSTRUME'].split('_')[0]
            proc_type = oimg[0].header['proctype']
        except KeyError:
            try:
                date = oimg[1].header['Date-obs']
                inst = oimg[1].header['INSTRUME'].split('_')[0]
                proc_type = oimg[1].header['proctype']
            except KeyError:
                print(f)
                raise Exception('Something Broke!')

    inst = inst.upper()
    date = date.split('T')[0] # we want the date not the time
    year = date.split('-')[0]
    month = date.split('-')[1]
    month = calendar.month_abbr[int(month)]

    target_dir = './{}{}{}/{}'.format(year, month, inst, proc_type.lower())

    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)

    # work out the relative path
    relpath = os.path.relpath(f, target_dir)
#    try:
#        os.symlink('{}'.format(relpath), '{}/{}'.format(target_dir, f))
#    except FileExistsError:
#        pass

    # move the fits image
    try:
        shutil.move(f, target_dir)
    except:
        pass

    # move the associated png if it exists
    try:
        shutil.move(f.replace('.fits.fz', '.png'), target_dir)
    except FileNotFoundError:
        pass
