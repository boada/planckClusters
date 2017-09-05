from astropy.io import fits
from glob import glob
import os

''' This file looks into the vo_proc directory and links everything into a
new directoy which is sorted by the object. We don't care when the
observations were taken but we care about what type of observation it is and
what istrument was used. All of the swarping and data processing is going to
happen in the proc directory.

'''

# use sorted files to keep the masks behind
files = sorted(glob('./vo_proc/201*/*/*.fz', recursive=True))

for i, f in enumerate(files):
    # keep the stacked/resampled files seperate
    # do the step first as it's the easiest to exclude files
    path_parts = f.split('/')
    if 'resampled' in f:
        file_type = 'resampled'
    elif 'stacked' in f:
        file_type = 'stacked'
        continue
    else:
        continue

    # don't take the original files
    if 'orig' in f:
        continue

    # open the image and start reading the properties
    with fits.open(f) as oimg:
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

        if 'mask' in obj.lower():
            continue
        elif 'exposure' in obj.lower():
            continue
        elif 'psz' not in obj.lower():
            continue

        print(obj)

    # deal with mosaic3
    if 'mosaic3' in inst.lower():
        inst = 'mosaic3'

    # for linking into the proc dir
    target_dir = './{}/{}'.format('proc2', obj)
    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)

    # work out the relative path
    relpath = os.path.relpath(f, target_dir)

    try:
        # image file
        os.symlink('{}'.format(relpath),
            '{}/{}'.format(target_dir, path_parts[-1]))
    except FileExistsError:
        pass

    # link the dqmask into the new directory
    if 'tu' in path_parts[-1].lower():
        # link the file name index +1
        new_file = files[i + 1]
        fname = new_file.split('/')[-1]
        #fname = path_parts[-1]
        #findex = fname.split('tu')[-1].rstrip('fits.fz')
        #findex_new = int(findex) + 1
        #fname = 'tu{}.fits.fz'.format(findex_new)
    elif 'k4' in path_parts[-1].lower():
        fname = path_parts[-1]
        if 'opi' in fname:
            fname = fname.replace('opi', 'opd')
        else:
            print('something is wrong!!! inner!!!')
    else:
        print('something is wrong!!! outer!!!')

    # update the file name
    path_parts[-1] = fname
    f = './{}'.format('/'.join(path_parts))
    # work out the relative path
    relpath = os.path.relpath(f, target_dir)

    try:
        # dqmask
        os.symlink('{}'.format(relpath),
            '{}/{}'.format(target_dir, fname))
    except FileExistsError:
        pass
