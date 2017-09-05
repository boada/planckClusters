import os
import sys
from astropy.io import ascii
import subprocess
import shlex

def webpage(tilename):

    bpzpath = os.environ['BPZPATH']

    if not os.path.isfile('{}_bpz.cat'.format(tilename)):
        cmd = 'python {}/bpzfinalize.py {}'.format(bpzpath, tilename)

        p = subprocess.Popen(shlex.split(cmd))
        p.wait()

    bpz = ascii.read('{}_bpz.cat'.format(tilename))
    mask = bpz['zspec'] != -99.
    bpz = bpz[mask]

    ids = ','.join([str(i) for i in bpz['id']])

    cmd = 'python {}/plots/webpage.py {} {}'.format(bpzpath, tilename, ids)
    print(cmd)

    p = subprocess.Popen(shlex.split(cmd))
    p.wait()


if __name__ == "__main__":
    webpage(sys.argv[1])

