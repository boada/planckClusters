from urllib import parse, request
import sys
import os
sys.path.append(f'{os.environ["HOME"]}/Projects/planckClusters/catalogs')
from load_catalogs import load_PSZcatalog
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
    select = '''SELECT p.objid,
                p.ra,
                p.raErr AS ra_err,
                p.dec,
                p.decErr AS dec_err,
                p.fiberMag_u,
                p.fiberMagErr_u,
                p.petroRad_u,
                p.petroRadErr_u,
                p.fiberMag_g,
                p.fiberMagErr_g,
                p.petroRad_g,
                p.petroRadErr_g,
                p.fiberMag_r,
                p.fiberMagErr_r,
                p.petroRad_r,
                p.petroRadErr_r,
                p.fiberMag_i,
                p.fiberMagErr_i,
                p.petroRad_i,
                p.petroRadErr_i,
                p.fiberMag_z,
                p.fiberMagErr_z,
                p.petroRad_z,
                p.petroRadErr_z,
                pz.z AS photoz,
                pz.zerr AS photoz_err,
                so.ra AS spec_ra,
                so.dec AS spec_dec,
                so.z AS specz,
                so.zerr AS specz_err,
                p.type,
                so.class
    '''
    FROM = '''FROM photoprimary AS p
    '''
    join = '''
       JOIN dbo.Fgetnearbyobjeq(''' + str(ra) + ',' + str(dec) + ''',30) AS N
         ON p.objid = N.objid
       LEFT JOIN photoz AS Pz
              ON p.objid = Pz.objid
       LEFT JOIN specobj AS SO
              ON p.objid = SO.bestobjid
              AND zwarning = 0
              AND SO.z > 0.01
    '''
    where = '''WHERE  p.i < 23
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
data = load_PSZcatalog()

for i, row in data.iterrows():
    n = row.NAME.replace(' ', '_')

    ra = float(row.RA)
    dec = float(row.DEC)

    work(ra, dec, f'./SDSS/{n}_SDSS_catalog.csv')
    sleep(1)
