from urllib import parse, request
import sys
import os
from astLib import astCoords
import numpy as np
from time import sleep

def filtercomment(sql):
    "Get rid of comments starting with --"
    fsql = ''
    for line in sql.split('\n'):
        fsql += line.split('--')[0] + ' ' + os.linesep
    return fsql


def query(sql):
    '''Run query and return file object'''

    url = 'http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
    fmt = 'csv'
    fsql = filtercomment(sql)
    params = parse.urlencode({'cmd': fsql, 'format': fmt})
    return request.urlopen(url + '?%s' % params)

def work(ra, dec, outfile):
    select = '''SELECT G.objid,
                G.ra,
                G.dec,
                G.raErr as ra_err,
                G.decErr as dec_err,
                G.cModelMag_u as u,
                G.cModelMagErr_u as u_err,
                G.cModelMag_g as g,
                G.cModelMagErr_g as g_err,
                G.cModelMag_r as r,
                G.cModelMagErr_r as r_err,
                G.cModelMag_i as i,
                G.cModelMagErr_i as i_err,
                G.cModelMag_z as z,
                G.cModelMagErr_z as z_err,
                G.flags,
                Pz.z AS Photoz,
                Pz.zerr AS Photoz_err,
                SO.specobjid,
                SO.ra AS Spec_ra,
                SO.dec AS Spec_dec,
                SO.z AS Specz,
                SO.zerr AS Specz_err
    '''
    FROM = '''FROM   galaxy AS G
    '''
    join = '''
       JOIN dbo.Fgetnearbyobjeq(''' + str(ra) + ',' + str(dec) + ''',30) AS GN
         ON G.objid = GN.objid
       LEFT JOIN photoz AS Pz
              ON G.objid = Pz.objid
       LEFT JOIN specobj AS SO
              ON G.objid = SO.bestobjid
              AND SO.zwarning = 0
    '''
    where = '''WHERE  G.i < 23
       AND clean = 1
       AND ( calibstatus_r & 1 ) != 0
    '''

    sql = select + FROM + join + where

    result = query(sql)

    line = result.readline()
    if line.startswith("ERROR".encode()):
        print('ERROR')
        ofp = sys.stderr
    else:
        with open(outfile, 'wt') as ofp:
            cnt = 0
            while line:
                ofp.write(line.rstrip().decode() + os.linesep)
                line = result.readline()
                cnt += 1
        if cnt < 3:
            os.remove(outfile)


# get file data
data = np.genfromtxt('../catalogs/PSZ2_unconfirmed_catalog - Master.csv',
           delimiter=',', names=True, dtype=None)

for i, (ra, dec,
        name) in enumerate(zip(data['RA'], data['DEC'], data['Name'])):
    print(data['Name'][i])

    ra = astCoords.hms2decimal(ra, ':')
    dec = astCoords.dms2decimal(dec, ':')
    work(ra, dec, './SDSS/%s_SDSS_catalog.csv' % (name.decode()))
    sleep(1)
