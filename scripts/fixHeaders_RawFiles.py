from astropy.io import fits
from glob import glob

files = sorted(glob('./raw/201*/*/*.fz', recursive=True))

for f in files:
    print(f)
    with fits.open(f,) as oimg:
        changed = False
        try:
            obj = oimg[0].header['object']
            if 'flat' in obj.lower():
                oimg[0].header['object'] = 'flat'
                changed = True
            elif 'test' in obj.lower():
                oimg[0].header['object'] = 'test'
                changed = True
            elif 'focus' in obj.lower():
                oimg[0].header['object'] = 'focus'
                changed = True
            elif 'dark' in obj.lower():
                oimg[0].header['object'] = 'dark'
                changed = True
            # mosaic specific stuff
            elif 'zero' in obj.lower():
                oimg[0].header['object'] = 'zero'
                changed = True
            elif 'standard' in obj.lower():
                oimg[0].header['object'] = 'standard'
                changed = True
            else:
                continue
        except IndexError:
            print('{} only has one extension'.format(f))
            continue
        except KeyError:
            print("{} isn't a image file?".format(f))
            continue

        if changed:
            #oimg.flush(output_verify='ignore')
            oimg.writeto(f, clobber=True, output_verify='ignore')
