from glob import glob
import os
import time
import sys

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


if __name__ == "__main__":
    try:
        if 'rebuild' in sys.argv[1].lower():
            rebuild = True
    except IndexError:
        rebuild = False
    dirs = glob('./PS*', recursive=True)
    cmds = mkCommands(dirs)
    mkTracker(cmds)
    if rebuild:
        rebuildTracker(cmds)
    else:
        doWork()
