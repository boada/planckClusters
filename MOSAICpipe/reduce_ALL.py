import os
from glob import glob

''' This file links the MOSAIC pipeline into each folder and then does the
complete reduction on things. It still needs to have the individual association
files created before hand, but it does everything else.

I've updated it to also to the newfirm linking and reduction. You specify which
instrument you want to use as a command line argument. 'mosaic' or 'newfirm'

'''

script_dir = '/home/boada/Projects/planckClusters/MOSAICpipe'

def main():
    dirs = [dirs for _, dirs, _ in os.walk('./')][0] # only want top level
    cwd = os.getcwd()
    for d in dirs:
        print(d)
        os.chdir(cwd)
        if 'PSZ' not in d:
            continue
        target_dir = './{}'.format(d)

        if not os.path.isdir(target_dir):
            continue

        relpath = os.path.relpath('{}'.format(script_dir), target_dir)
        print(relpath)
        print(target_dir)
        try:
            os.symlink('{}/combcat_PROJECTED.py'.format(script_dir),
                    '{}/combcat_PROJECTED.py'.format(target_dir))
        except FileExistsError:
            pass

        # now do the pipeline
        os.chdir(target_dir)

        assocFile = glob('*.assoc')[0]
        print(os.getcwd())
        # build the command
        cmd = 'python3 combcat_PROJECTED.py {} ./ ./'.format(assocFile)
        cmd += ' --bpz --dust --newfirm'

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
    main()
