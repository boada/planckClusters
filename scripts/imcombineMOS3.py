from glob import glob
import sys
import tempfile
import os
import pyraf.iraf
from itertools import repeat

def imcombine(images, output_path, flow, fhigh, **kwargs):
    """ Mim-max average combine all the images in the set.
    The method invokes the IRAF task 'imcombine' on the FITS images contained
    in 'images', an iterable. The type of combining operation performed on the
    pixels (the 'combine' parameter when 'imcombine' is invoked) is 'average',
    while the rejection operation (the 'reject' parameter) is 'minmax'. That is
    the reason why 'flow' and 'fhigh', which determine the fraction of low and
    high pixels, respectively, that are to be rejected, are required as
    arguments. If existing, the output image is silently overwritten.
    The WCS in the image is used to derive the offsets.
    """

    if not len(images):
        raise ValueError("no FITS images given")

    if not 0 <= flow <= 1 or not 0 <= fhigh <= 1:
        raise ValueError("both 'nlow' and 'fhigh' must be in the range [0,1]")

    # Number of low and high pixels are to reject; round to nearest integer
    nlow = int(round(len(images) * flow))
    nhigh = int(round(len(images) * fhigh))

    # If only two images are combined, we always discard the highest pixel
    # (since some stars could be discernible in the case of flat-fields,
    # for example). When there are three or more images we also discard the
    # lowest and highest pixels, unless the specified fraction is zero: in
    # that case, it is understood as to mean "do not discard anything".
    if len(images) == 2 and not nhigh and fhigh:
        nhigh = 1
    elif len(images) >= 3:
        if not nlow and flow:
            nlow = 1
        if not nhigh and fhigh:
            nhigh = 1

    # There is what seems to be a bug in PyRAF that makes it impossible to
    # imcombine more than approximately forty images, so we save the list
    # of input images to a temporary file, use it as the input to IRAF's
    # imcombine and delete it before leaving the method.

    input_fd, input_path = tempfile.mkstemp(suffix='.lst', text=True)

    try:
        for path in images:
            os.write(input_fd, '{0}\n'.format(path).encode())
        os.close(input_fd)

        if os.path.exists(output_path):
            os.unlink(output_path)

        pyraf.iraf.imcombine('@' + input_path,
                             output_path,
                             combine='sum',
                             reject='minmax',
                             nlow=nlow,
                             nhigh=nhigh,
                             offsets='wcs',
                             **kwargs)
    finally:
        os.unlink(input_path)

def fix_header(image, mask=False):
    pyraf.iraf.hedit(image,
                     fields='PV[12]_*',
                     delete='yes',
                     verify='no')
    if mask:
        pyraf.iraf.hedit(image, fields='PRODTYPE', value='dqmask',
                         add='yes', verify='no')


if __name__ == "__main__":

    if os.path.isdir(sys.argv[1]):
        files = glob('{}/*.fits'.format(sys.argv[1]))
        for f in files:
            print(f)
            input_paths = ['{}[{}]'.format(f, i + 1) for i, f in
                           enumerate(repeat(f, 4))]
            output_path = f.rstrip('.fits') + '_multi.fits'
            if 'opi' in f:
                imcombine(input_paths, output_path, 0.0, 0.0, blank=0)
                fix_header(output_path)
            elif 'opd' in f:
                imcombine(input_paths, output_path, 0.0, 0.0, blank=1,
                          outtype='integer')
                fix_header(output_path, mask=True)
    else:
        input_paths = sys.argv[1:-1]
        output_path = sys.argv[-1]
        imcombine(input_paths, output_path, 0.0, 0.0)
        fix_header(output_path)
