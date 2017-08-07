import pandas as pd

base_string = ('https://ned.ipac.caltech.edu/cgi-bin/objsearch?in_csys='
               'Equatorial&in_equinox=J2000.0&lon={:.6f}d&lat={:.6f}d'
               '&radius=5&hconst=73&omegam=0.27&omegav=0.73&corr_z=1'
               '&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&in_'
               'objtypes1=Galaxies&in_objtypes1=GPairs&in_objtypes1=GTriples&'
               'in_objtypes1=GGroups&in_objtypes1=GClusters&search_type='
               'Near+Position+Search&nmp_op=ANY&out_csys=Equatorial&out_'
               'equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text'
               '&zv_breaker=30000.0&list_limit=5&img_stamp=YES')

sheet = pd.read_csv('./catalogs/PSZ2_unconfirmed_catalog - Master.csv')

for index, row in sheet.iterrows():
    sheet['NED'][index] = base_string.format(row['RA'], row['DEC'])

sheet.to_csv('updated.csv')
