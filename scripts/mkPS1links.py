import os
import sys
sys.path.append(f'{os.environ["HOME"]}/Projects/planckClusters/catalogs')
from load_catalogs import load_PSZcatalog

base_string = ('http://ps1images.stsci.edu/cgi-bin/ps1cutouts?'
               'pos={}%2C{}&'
               'filter=color&filter=g&filter=r&filter=i&filetypes=stack&'
               'auxiliary=data&size=1200&output_size=0&verbose=0&'
               'autoscale=99.500000&catlist=')

data = load_PSZcatalog(unconf=True)
data['PanSTARRS'] = ''

for index, row in data.iterrows():
    if row.DEC < -31:
        continue
    else:
        data['PanSTARRS'][index] = base_string.format(row.RA, row.DEC)
data.to_csv('updated_withPS1.csv', index=False)








