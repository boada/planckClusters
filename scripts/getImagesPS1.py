import numpy
from time import sleep
import subprocess
import shlex
from astLib.astCoords import dms2decimal
def work(ra, dec, name):

    cmd = 'panstamps --width=10 '
    cmd += '--downloadFolder=./PS1/{} stack {} {}'.format(name, ra, dec)
    print(cmd)
    p = subprocess.Popen(shlex.split(cmd))

    p.wait()


data = numpy.genfromtxt('../catalogs/PSZ2_unconfirmed_catalog - Master.csv',
                     delimiter=',', names=True, dtype=None)

for i, (ra, dec,
        name) in enumerate(zip(data['RA_SEX'], data['DEC_SEX'], data['Name'])):
    print(data['Name'][i])

    dec_decimal = dms2decimal(dec.decode(), ':')

    if dec_decimal < -31:
        continue

    #name = ''.join(e for e in name.decode() if e.isalnum())
    work(ra.decode(), dec.decode(), '{}'.format(name.decode()))
    sleep(1)
