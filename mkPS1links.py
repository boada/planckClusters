import pandas as pd

base_string_pos = ('http://ps1images.stsci.edu/cgi-bin/ps1cutouts?'
               'pos={}%3A{}%3A{}{}%2B{}%3A{}%3A{}&'
               'filter=color&filter=g&filter=r&filter=i&filetypes=stack&'
               'auxiliary=data&size=1200&output_size=0&verbose=0&'
               'autoscale=99.500000&catlist=')

base_string_neg = ('http://ps1images.stsci.edu/cgi-bin/ps1cutouts?'
               'pos={}%3A{}%3A{}+{}{}%3A{}%3A{}&'
               'filter=color&filter=g&filter=r&filter=i&filetypes=stack&'
               'auxiliary=data&size=1200&output_size=0&verbose=0&'
               'autoscale=99.500000&catlist=')

sheet = pd.read_csv('./catalogs/PSZ2_unconfirmed_catalog - Master.csv')

for index, row in sheet.iterrows():
    print(row['RA_SEX'], row['DEC_SEX'])
    try:
        ra_parts = row['RA_SEX'].split(':')
        dec_parts = row['DEC_SEX'].split(':')
        if int(dec_parts[0]) > 0:
            sheet['PanSTARRS'][index] = base_string_pos.format(ra_parts[0],
                                                    ra_parts[1],
                                                    ra_parts[2],
                                                    dec_parts[0][0], # the +/-
                                                    dec_parts[0][1:],
                                                    dec_parts[1],
                                                    dec_parts[2])
        if int(dec_parts[0]) < 0:
            sheet['PanSTARRS'][index] = base_string_neg.format(ra_parts[0],
                                                    ra_parts[1],
                                                    ra_parts[2],
                                                    dec_parts[0][0], # the +/-
                                                    dec_parts[0][1:],
                                                    dec_parts[1],
                                                    dec_parts[2])
    except:
        pass
sheet.to_csv('updated_withPS1.csv')








