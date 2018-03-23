from astropy.io import fits
from glob import glob
from numpy import genfromtxt, where

''' This makes sure that the names of the objects are consistent. The files
downloaded from the VO sometimes have typos in the names or inconsistent
spaces. This little script modifies the header keywords to make sure that
everything is the same throughout. Most of the changes are at the bottom and
most of them come by looking at the directory names in the proc directory. If
you see something funny in there then you can add a rule to this script which
will fix it.

Really this should only have to be run one time and then it will fix everything.

'''


# use sorted files to keep the masks behind
files = sorted(glob('./vo_proc/201*/*/*.fz', recursive=True))

data = genfromtxt('/home/boada/Projects/planckClusters/'
                'catalogs/PSZ2_unconfirmed_catalog - current.csv',
                  delimiter=',',
                  names=True, dtype=None)
for i, f in enumerate(files):
    # keep the stacked/resampled files seperate
    # do the step first as it's the easiest to exclude files
    path_parts = f.split('/')
    if 'resampled' in f:
        file_type = 'resampled'
    elif 'stacked' in f:
        file_type = 'stacked'
    elif 'opd' in f:
        continue
    else:
        continue

    obj_old = ''
    # open the image and start reading the properties
    with fits.open(f, mode='update') as oimg:
        try:
            obj = oimg[0].header['object']
            inst = oimg[0].header['instrume'].split('_')[0]
        except KeyError:
            try:
                obj = oimg[1].header['object']
                inst = oimg[1].header['instrume'].split('_')[0]
            except KeyError:
                print('still have problem')
                obj = 'unknown'

        # fix the filter info if there is a problems
        try:
            for hdr, _ in enumerate(oimg):
                if oimg[hdr].header['FILTER'] == 'z SDSS k1020':
                    print('wrong filter!')
                    oimg[hdr].header['FILTER'] = 'z SDSS c6020'
                    oimg[hdr].header['FILTID'] = 'c6020'
        except KeyError:
            pass

        if 'mask' in obj.lower():
            continue
        elif 'exposure' in obj.lower():
            continue

        obj_old = obj

        # here is where we fix the header object keyword
        # lets fix the specific problems first
        if obj == 'RU 149F' or obj == 'Ru149':
            obj = 'Ru149F'
        elif obj == 'SA 92-282':
            obj = '92_282'
        elif obj == 'SA 95-193':
            obj = '95_193'
        elif obj == 'Hilt 1089':
            obj = 'Hilt1089'
        elif obj == 'BD+25 1981':
            obj = 'BD+25o1981'
        elif obj == 'PSZ1_G153.56+36.32':
            obj = 'PSZ1_G153.56+36.23'
        elif obj == 'PSZ2_G093.71-30.91':
            obj = 'PSZ2_G093.71-30.90'
        elif obj == 'PSZ2_203.32+08.91':
            obj = 'PSZ2_G203.32+08.91'
        elif obj == 'OBJECT OBSERVATION(s)':
            obj = 'PSZ2_G210.71+63.08'

        # one set of images has two clusters in it... Sometimes they were split
        # and sometimes they weren't. This combines them all into a single
        # thing.
        elif obj == 'PSZ2_G153.56+PSZ2_G153.68' or obj == 'PSZ2_G153.68+36.96':
            obj = 'PSZ2_G153.56+36.82'

        # now we fix the extra spaces
        elif ' ' in obj:
            obj = obj.replace(' ', '_')
        # commas instead of periods
        elif ',' in obj:
            obj = obj.replace(',', '.')
        # two many periods
        elif '..' in obj:
            obj = obj.replace('..', '.')
        # typos
        elif 'PSZ_' in obj:
            obj = obj.replace('PSZ_', 'PSZ2_')

        # now we need to change all of the PSZ1 names to PSZ2 names. Going to
        # use the catalog we loaded at the top of this block
        elif 'PSZ1' in obj:
            try:
                indx = where(data['PSZ1_Name'] == obj.encode())[0][0]
                obj = data['Name'][indx].decode() # this is the PSZ2 name
            except IndexError:
                pass

        if obj_old != obj:
            print('updated: {} --> {}'.format(obj_old, obj))
            oimg[0].header['object'] = obj
            oimg[1].header['object'] = obj
