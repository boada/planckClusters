import os
import numpy as np
from time import sleep

def work(ra, dec, target):
    sql = '''SELECT o.objid,
              o.ramean,
              o.decmean,
              m.gmeanpsfmag,
              m.gmeanpsfmagerr,
              m.rmeanpsfmag,
              m.rmeanpsfmagerr,
              m.imeanpsfmag,
              m.imeanpsfmagerr,
              m.imeankronmag,
              m.zmeanpsfmag,
              m.zmeanpsfmagerr,
              m.ymeanpsfmag,
              m.ymeanpsfmagerr,
              m.gmeanpsfmagnpt,
              m.rmeanpsfmagnpt,
              m.imeanpsfmagnpt,
              m.zmeanpsfmagnpt,
              m.gflags,
              m.gqfperfect,
              m.rflags,
              m.rqfperfect,
              m.iflags,
              m.iqfperfect,
              m.zflags,
              m.zqfperfect
            INTO   MyDB.{}
            FROM   Fgetnearbyobjeq({}, {}, 10) AS nb
                INNER JOIN objectthin AS o
                        ON o.objid = nb.objid
                            AND o.ndetections > 1
                INNER JOIN meanobject AS m
                        ON o.objid = m.objid
                            AND o.uniquepspsobid = m.uniquepspsobid
            WHERE  m.gqfperfect >= 0.85
            AND m.rqfperfect >= 0.85
            AND m.iqfperfect >= 0.85
            AND m.zqfperfect >= 0.85
            AND m.imeanpsfmag - m.imeankronmag > 0.5
    '''.format(target, ra, dec)

    # build the command
    cmd = 'java -jar casjobs.jar run "{}"'.format(sql)

    #os.system(cmd)
    print(cmd)

    cmd = 'java -jar casjobs.jar extract -b '
    cmd += '{} -F -d PS1/'.format(target)

    #os.system(cmd)
    print(cmd)

    cmd = 'java -jar casjobs.jar execute -t "MyDB/1" -n '
    cmd += '"drop query" "drop table {}"'.format(target)

    #os.system(cmd)
    print(cmd)

# get file data
data = np.genfromtxt('../catalogs/PSZ2_unconfirmed_catalog - proc2.csv',
           delimiter=',', names=True, dtype=None)

for i, (ra, dec,
        name) in enumerate(zip(data['RA'], data['DEC'], data['Name'])):
    print(data['Name'][i])

    if dec < -31:
        continue

    name = ''.join(e for e in name.decode() if e.isalnum())
    work(ra, dec, '{}'.format(name))
    sleep(1)
    break
