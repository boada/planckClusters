import urllib
import sys
import os
import string
from astLib import astCoords
import numpy as np
import glob
import time


def filtercomment(sql):
    "Get rid of comments starting with --"
    fsql = ''
    for line in sql.split('\n'):
        fsql += line.split('--')[0] + ' ' + os.linesep
    return fsql


def query(sql):
    url = 'http://skyserver.sdss3.org/public/en/tools/search/x_sql.aspx'
    fmt = 'csv'
    "Run query and return file object"
    fsql = filtercomment(sql)
    params = urllib.urlencode({'cmd': fsql, 'format': fmt})
    return urllib.urlopen(url + '?%s' % params)


def work(ra, dec, outfile):
    select = '''SELECT TOP 1000 G.objid,
                G.ra,
                G.dec,
                G.u,
                G.err_u,
                G.g,
                G.err_g,
                G.r,
                G.err_r,
                G.i,
                G.err_i,
                G.z,
                G.err_z,
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
       JOIN dbo.Fgetnearbyobjeq(''' + str(ra) + ',' + str(dec) + ''',4) AS GN
         ON G.objid = GN.objid
       LEFT JOIN photoz AS Pz
              ON G.objid = Pz.objid
       LEFT JOIN specobj AS SO
              ON G.objid = SO.bestobjid
    '''
    where = '''WHERE  G.r < 24
       AND clean = 1
       AND ( calibstatus_r & 1 ) != 0
    '''

    sql = select + FROM + join + where

    result = query(sql)

    line = result.readline()
    if line.startswith("ERROR"):
        ofp = sys.stderr
    else:
        ofp = open(outfile, 'wt')
        while line:
            ofp.write(string.rstrip(line) + os.linesep)
            line = result.readline()
        ofp.close()

# get file data
data =\
np.genfromtxt('../../PSZ2_unconfirmed_catalog_4NOAO_2016A_newSwift-PSZ2_unconfirmed_catalog_4NOAO_.csv',
           delimiter=',', names=True, dtype=None)

for i, (ra, dec,
        name) in enumerate(zip(data['RA'], data['Dec'], data['Name'])):
    if not data['SDSS_Footprint'][i] == 'TRUE':
        continue
    print(data['Name'][i])
    if not os.path.isdir(data['Name'][i]):
        os.mkdir(data['Name'][i])

    ra = astCoords.hms2decimal(ra, ':')
    dec = astCoords.dms2decimal(dec, ':')
    work(ra, dec, './%s/%s_SDSS_catalog.txt' % (name, name))
