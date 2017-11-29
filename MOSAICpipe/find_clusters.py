from glob import glob
import os
import time


def mkCommands(dirs):
    cmds = []
    for d in dirs:
        parts = d.split('/') # get the path parts

        cluster = parts[1]

        # start building the command
        cmd = 'python3 find_position.py'

        # build the rest
        cmd += ' {} --path {} --dx 1000 --dy 1000'.format(cluster, d)

        cmd += ' --pixelscale 0.25'

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

    with open('tracker', 'w') as tracker:
        for i, (done, cmd) in enumerate(zip(status, cmds)):
            if done == 'True':
                continue
            if done == 'False':
                print(cmd)
                os.system(cmd)
                status[i] = 'True'
            print('sleeping for 5 seconds. Press ctrl-c to exit')
            try:
                time.sleep(5)
            except KeyboardInterrupt:
                print('updating tracker...')
                for done, cmd in zip(status, cmds):
                    tracker.write('{} {}\n'.format(done, cmd))
                return


dirs = glob('./PS*', recursive=True)
cmds = mkCommands(dirs)
mkTracker(cmds)
doWork()
