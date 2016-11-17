import montage_wrapper as montage
import os
import multiprocessing
import numpy as np
from glob import glob
import traceback, functools

def error(msg, *args):
    multiprocessing.log_to_stderr()
    return multiprocessing.get_logger().error(msg, *args)

def trace_unhandled_exceptions(func):
    @functools.wraps(func)
    def wrapped_func(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except Exception as e:
            error(traceback.format_exc())
            raise

    return wrapped_func

class AsyncFactory:
    def __init__(self, func, cb_func):
        self.func = func
        self.cb_func = cb_func
        self.pool = multiprocessing.Pool(maxtasksperchild=1,
                                         processes=multiprocessing.cpu_count())

    def call(self,*args, **kwargs):
        self.pool.apply_async(self.func, args, kwargs, self.cb_func)

    def wait(self):
        self.pool.close()
        self.pool.join()

@trace_unhandled_exceptions
def worker(f):
    print("PID: %d \t Cluster: %s" % (os.getpid(), f))
    os.chdir(f)
    for band in 'ugriz':
        if not os.path.isdir(band):
            os.mkdir(band)
            images = glob('*_{}_*.fits'.format(band))
            for i in images:
                os.rename('./{}'.format(str(i)), './{}/{}'.format(band,
                                                                  str(i)))

        # make the mosaics
        try:
            montage.mosaic(band, '{}_mosaic'.format(band), background_match=True,
                       combine='median', work_dir='./{}_tmp'.format(band))
        except OSError as e:
        # montage throws an error if any directory already exists. Check all
        # directories to make sure we haven't missed anything.
            error(traceback.format_exc())
            continue

    return f

def cb_func(f):
    print("PID: %d \t Value: %s completed" % (os.getpid(), f))

if __name__ == "__main__":

    async_worker = AsyncFactory(worker, cb_func)
    # get file data

    files = glob('PSZ2*')

    print(len(files))
    for i, f in enumerate(files):
        async_worker.call(f)

    async_worker.wait()

