import os
from glob import glob
import sys

''' This file links the MOSAIC pipeline into each folder and then does the
complete reduction on things. It still needs to have the individual association
files created before hand, but it does everything else.

I've updated it to also to the newfirm linking and reduction. You specify which
instrument you want to use as a command line argument. 'mosaic' or 'newfirm'

Basically just specify mosaic or mosaic3 and this will do the rest. You will
need to do both mosaic and mosaic3 but it doesn't matter which you do first.

'''

script_dir = '/home/boada/Projects/planckClusters/MOSAICpipe'

def main(inst='mosaic'):
    dirs = [dirs for _, dirs, _ in os.walk('./')][0] # only want top level
    cwd = os.getcwd()
    for d in dirs:
        print(d)
        os.chdir(cwd)
        if 'PSZ' not in d:
            continue
        target_dir = './{}/{}/resampled'.format(d, inst)

        if not os.path.isdir(target_dir):
            continue

        relpath = os.path.relpath('{}'.format(script_dir), target_dir)
        print(relpath)
        print(target_dir)
        try:
            if inst == 'mosaic':
                os.symlink('{}/merge_catalogs_MOS1.py'.format(script_dir),
                        '{}/merge_catalogs_MOS1.py'.format(target_dir))
            elif inst == 'mosaic3':
                os.symlink('{}/merge_catalogs_MOS1.py'.format(script_dir),
                        '{}/merge_catalogs_MOS1.py'.format(target_dir))
            else:
                print('instrument not understood')
                return
        except FileExistsError:
            pass

        # now do the pipeline
        os.chdir(target_dir)

        assocFile = glob('*.assoc')[0]
        print(os.getcwd())
        # build the command
        if inst == 'mosaic':
            cmd = 'python3 merge_catalogs_MOS1.py {} ./'.format(
                assocFile.rstrip('.assoc'))
        elif inst == 'mosaic3':
            cmd = 'python3 merge_catalogs_MOS1.py {} ./'.format(
                assocFile.rstrip('.assoc'))
        else:
            return

        print(cmd)
        os.system(cmd)

    # clean up all of the intermediate data products
    cmds = ["find . -path '*/.diagnostics/*' -delete",
            "find . -type d -name '.diagnostics' -empty -delete",
            "find . -type f -name 'registration_*' -delete",
            "find . -type f -name '*ldac*' -delete",
            "find . -type f -name 'diagnostics.html' -delete",
            "find . -type f -name '*.lst' -delete",
            "find . -type f -name '*.xml' -delete",
            "find . -type f -name 'GAIA.cat' -delete",
            "find . -type f -name 'best_astrometry.dat' -delete"]
    for cmd in cmds:
        os.system(cmd)


if __name__ == "__main__":
    main(sys.argv[1])
