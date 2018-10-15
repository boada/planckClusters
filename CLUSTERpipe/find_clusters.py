from glob import glob
import os
import time
from numpy import genfromtxt, where

def mkCommands(dirs, radec=False):
    cmds = []
    for d in dirs:
        parts = d.split('/') # get the path parts

        cluster = parts[1]
        # start building the command
        cmd = 'python3 find_position.py'

        # build the rest
        cmd += ' {} --path {} '.format(cluster, d)
        cmd += ' --pixelscale 0.25'

        if radec:
            # patch for a silly situation
            if cluster == 'PSZ2_G128.15-24.17':
                cluster = 'PSZ2_G128.15-24.71'
            cat_dir = '/home/boada/Projects/planckClusters/catalogs'
            data = genfromtxt(cat_dir + '/PSZ2_unconfirmed_catalog - current.csv',
                              delimiter=',', names=True, dtype=None)

            # find the ra and dec position
            # find the index first -- the zeros just remove the arrays
            indx = where(data['Name'] == cluster.encode())[0][0]

            cmd += ' --RA {}'.format(data['RA_SEX'][indx].decode())
            cmd += ' --DEC {}'.format(data['DEC_SEX'][indx].decode())

        cmds.append(cmd)

    return cmds

def mkTracker(cmds, clobber=False):
    if not os.path.isfile('tracker') or clobber:
        with open('tracker', 'w') as tracker:
            for cmd in cmds:
                tracker.write('False {}\n'.format(cmd))

        return
    else:
        print('Tracker already exists. Call with clobber=True to overwrite')
        return

def rebuildTracker(cmds, clobber=True):
    imgs = glob('PSZ*/**/*A.png', recursive=True)
    fields = [i.split('/')[0] for i in imgs] # the fields already done
    if not os.path.isfile('tracker') or clobber:
        with open('tracker', 'w') as tracker:
            for cmd in cmds:
                status = False
                for i, f in enumerate(fields):
                    if f in cmd:
                        status = True
                        break
                if status:
                    tracker.write('True {}\n'.format(cmd))
                else:
                    tracker.write('False {}\n'.format(cmd))

        return

def doWork():
    print('doing work....')

    status = []
    cmds = []
    # read the whole file... it's not that big.
    with open('tracker', 'r') as tracker:
        for line in tracker:
            done = line.split()[0]
            cmd = line.split()[1:]
            # rebuild the command
            cmd = ' '.join(cmd)

            status.append(done)
            cmds.append(cmd)

    fields = len(cmds)
    with open('tracker', 'w') as tracker:
        for i, (done, cmd) in enumerate(zip(status, cmds)):
            if done == 'True':
                print()
                print()
                print('{}/{} Fields completed.'.format(i + 1, fields))
                print('Sleeping for 5 seconds. Press ctrl-c to exit')
                continue
            if done == 'False':
                print(cmd)
                os.system(cmd)
                status[i] = 'True'
            print()
            print()
            print('{}/{} Fields completed.'.format(i + 1, fields))
            print('Sleeping for 5 seconds. Press ctrl-c to exit')
            try:
                time.sleep(5)
            except KeyboardInterrupt:
                print('updating tracker...')
                for done, cmd in zip(status, cmds):
                    tracker.write('{} {}\n'.format(done, cmd))
                return

def cmdline():
    from optparse import OptionParser

    USAGE = 'usage:\t % prog [options] \n'
    USAGE += 'i.e.: %prog --rebuild or %prog --radec'

    parser = OptionParser(usage=USAGE)

    parser.add_option('--rebuild', dest='rebuild', default=False,
                      help='Rebuild Tracker')
    parser.add_option('--radec', dest='radec', default=False,
                      action='store_true',
                      help='Use RA/DEC for centering')

    options, args = parser.parse_args()

    return options, args


if __name__ == "__main__":

    opt, args = cmdline()

    # get the directory list
    dirs = glob('./PS*', recursive=True)

    cmds = mkCommands(dirs, radec=opt.radec)
    mkTracker(cmds)

    if opt.rebuild:
        rebuildTracker(cmds)
    else:
        doWork()
