#!/usr/bin/env python
""" Make a ds9 region file from a SExtractor catalogue."""
from __future__ import division, print_function, unicode_literals
try:
    unicode
except NameError:
    import pickle
    unicode = basestring = str
    xrange = range
else:
    import cPickle as pickle
import numpy as np

usage = """\
Usage: sex2DS9reg sextractor_catalogue_filename [DS9region_filename]
Create a DS9 region file from a SExtractor output catalogue. It must
contain either X_IMAGE and Y_IMAGE, or XWIN_IMAGE and YWIN_IMAGE. If
A, B and THETA columns are also present, then ellipse regions are
drawn, otherwise just points.
"""


def readtxt(fh,
            sep=None,
            usecols=None,
            comment='#',
            skip=0,
            arrays=True,
            names=None,
            readnames=False,
            converters=None,
            mintype=int):
    """ Reads columns from a text file into arrays, converting to int,
    float or str where appropriate.
    By default the column separator is whitespace. `rows` can be
    either an input filename or an iterable (e.g. a file object, list
    or iterator).
    Parameters
    ----------
    rows : filename or iterable object
        Input data.
    sep : str (default `None`)
        A string used to separate items on a row (also known as a
        delimiter). Default is None, which means whitespace.
    usecols : int or tuple of ints, optional
        Indices of columns to be read. By default all columns are read.
    comment : str (default `#`)
        Character marking the start of a comment.
    skip : int (default `0`)
        Number of rows to skip (not counting commented or blank lines)
        before reading data.
    arrays : bool (`True`)
        If True, all columns are converted to Numpy arrays.  If False,
        columns are returned as lists.
    names : str or sequence of str (default `None`)
        If `names` is given and `arrays` is True, the data are placed
        in a Numpy record array with field names given by `names`. Can
        also be a single string of comma-separated values.
    readnames : bool (`False`)
        If `readnames` is True the first line of the file is read to
        find the field names. This overrides the `names` keyword.
    converters : dict (`None`)
        Functions to apply to each entry of a column. Each (key,value)
        pair gives the column index (key) and the function to be
        applied to each entry in that column (value).
    Returns either structured array or lists.
    Examples
    --------
    >>> list_of_all_cols = readtxt('filename')
    >>> ninthcol, fifthcol = readtxt('filename', sep=',', usecols=(8,4)])
    >>> firstcol = readtxt('filename', comment='%', usecols=[0])
    >>> recarray = readtxt('filename', sep=',', usecols=(1,3), names='x,y'])
    """
    if mintype == float:
        typedict = {float: lambda x: str(x).strip()}
    elif mintype == int:
        typedict = {int: float, float: lambda x: str(x).strip()}
    else:
        raise ValueError('Unknown minimum type %s' % mintype)

    def convert(row, funcs):
        # convert each item in a row to int, float or str.
        for i, item in enumerate(row):
            while True:
                try:
                    row[i] = funcs[i](item)
                except ValueError:
                    # update the list of converters
                    try:
                        funcs[i] = typedict[funcs[i]]
                    except KeyError:
                        raise ValueError('Converter %s failed '
                                         'on %r' % (funcs[i], item))
                else:
                    break
        return row, funcs

    needclose = False
    if isinstance(fh, basestring):
        if fh.endswith('.gz'):
            import gzip
            fh = gzip.open(fh)
        else:
            fh = open(fh)
        needclose = True

    data = iter(fh)

    if comment is not None:
        len_comment = len(comment)

    if names and isinstance(names, basestring):
        names = [n.strip() for n in str(names).split(',')]
    elif names:
        names = map(str, names)

    skipped = 0
    out = []
    # main loop to read data
    for irow, row in enumerate(data):
        if comment is not None:
            row = row.split(comment)[0]
        row = row.lstrip()
        if not row: continue
        if skipped < skip:
            skipped += 1
            continue
        row = row.split(sep)
        if readnames:
            names = [str(r.strip()) for r in row]
            readnames = False
            continue
        if not out:
            # first row with data, so initialise converters
            funcs = [mintype] * len(row)
            if converters is not None:
                for i in converters:
                    funcs[i] = converters[i]
            if usecols is not None:
                funcs = [funcs[i] for i in usecols]
        if usecols is not None:
            try:
                row = [row[i] for i in usecols]
            except IndexError:
                raise IndexError('Columns indices: %s, but only %i entries in '
                                 'this row!' % (usecols, len(row)))
        try:
            row, funcs = convert(row, funcs)
        except IndexError:
            # Probably there are more items in this row than in
            # previous rows. This usually indicates a problem, so
            # raise an error.
            raise IndexError('Too many items on row %i: %s' % (irow + 1, row))

        if names:
            assert len(row) == len(names), '%i, %i, %s ' % (len(names),
                                                            irow + 1, row)
        out.append(row)

    if needclose:
        fh.close()

    # rows to columns, truncating to number of words on shortest line.
    if arrays:
        if names is not None:
            out = np.rec.fromrecords(out, names=names)
        else:
            out = [np.array(c) for c in zip(*out)]
    else:
        out = [list(c) for c in zip(*out)]

    if len(out) == 1 and names is None:
        return out[0]
    else:
        return out


def readsex(filename, catnum=None):
    """ Read a sextractor catalogue into a Numpy record array.
    Parameters
    ----------
    filename : str
      Sextractor output catalogue name
    catnum : int, optional
      If the Sextractor file is in LDAC_FITS format and contains more
      than one catalogue, this option specifies the catalogue number.
    Returns
    -------
    s : numpy record array
      Record array with field names the same as those in the
      sextractor catalogue.
    """
    fh = open(filename)
    # get the header
    row = fh.next()
    while not row.strip():
        row = fh.next()
    if row[8] == '=':
        fh.close()
        # assume a fits file
        try:
            import pyfits
        except ImportError:
            import astropy.io.fits as pyfits
        fh = pyfits.open(filename)
        if len(fh) > 3 and catnum is None:
            raise ValueError("specify catalogue number")
        elif catnum is not None:
            return pyfits.getdata(filename, catnum * 2).view(np.recarray)
        else:
            return pyfits.getdata(filename, 2).view(np.recarray)
    hd = []
    while row.startswith('#'):
        if row[1:].strip():
            hd.append(row)
        row = fh.next()
    fh.close()
    # get column numbers and names
    number, names = zip(*[row.split() for row in hd])[1:3]
    indcol = [int(c) - 1 for c in number]
    if len(names) - len(set(names)):
        dup = [n for n in set(names) if names.count(n) > 1]
        raise ValueError('fields with same names: %s' % dup)
    # read in the data
    return readtxt(filename, names=names, usecols=indcol)


def sex_to_DS9reg(filename,
                  s,
                  colour='green',
                  tag='all',
                  withtext=False,
                  use_WORLD=False):
    """Write a DS9 region file from SExtractor output.
    Parameters
    ----------
    filename : str
      Region file name.
    s : array
      The output of `readsex`.
    colour : str ('green')
      Region colour. One of {cyan blue magenta red green yellow white
      black}
    tag : str ('all')
      DS9 tag for all the regions
    with_text : bool (False)
      If True, then mark each region with either its magnitude (if
      given), otherwise its index in the input array `s`.
    use_WORLD : bool (False)
      If True, use WORLD coordinates (typically RA/Dec) instead of
      IMAGE coordinates.
    """

    names = set(s.dtype.names)
    regions = ['global font="helvetica 10 normal" select=1 highlite=1 '
               'edit=0 move=1 delete=1 include=1 fixed=0 source']
    fmt = 'ellipse(%s %s %s %s %s) # text={%s} color=%s tag={%s}'
    if not use_WORLD:
        regions.append('image')
        fields = ['X_IMAGE', 'Y_IMAGE']
        if not ('X_IMAGE' in names and 'Y_IMAGE' in names):
            fields = ['XWIN_IMAGE', 'YWIN_IMAGE']
            if not ('XWIN_IMAGE' in names and 'YWIN_IMAGE' in names):
                raise ValueError('require either X_IMAGE and Y_IMAGE '
                                 'or XWIN_IMAGE and YWIN_IMAGE')

        ellipse_vals = ['A_IMAGE', 'B_IMAGE', 'THETA_IMAGE']
        ellipsewin_vals = ['AWIN_IMAGE', 'BWIN_IMAGE', 'THETAWIN_IMAGE']
    else:
        regions.append('J2000')
        fields = ['X_WORLD', 'Y_WORLD']
        ellipse_vals = ['A_WORLD', 'B_WORLD', 'THETA_WORLD']

    if all((n in names) for n in ellipse_vals):
        fields = list(fields) + ellipse_vals
    elif all((n in names) for n in ellipsewin_vals):
        fields = list(fields) + ellipsewin_vals
    else:
        # we don't have any ellipticity info, just write points.
        fmt = 'point(%s %s) # point=circle text={%s} color=%s tag={%s}'

    for i, rec in enumerate(s):
        vals = [rec[f] for f in fields]
        if withtext:
            if 'MAG_AUTO' in names:
                text = '%i %.2f' % (i + 1, rec['MAG_AUTO'])
            else:
                text = i + 1
        else:
            text = ''
        if use_WORLD:
            vals[-1] = -vals[-1]
        vals.extend([text, colour, tag])
        regions.append(fmt % tuple(vals))

    fh = open(filename, 'w')
    fh.write('\n'.join(regions))
    fh.close()


def main(args):

    if len(args) not in (1, 2):
        print(usage)
        sys.exit()

    catname = args[0]
    if len(args) == 1:
        regname = catname.rsplit('.')[0] + '.reg'
    else:
        regname = args[1]

    print('Reading', catname)
    s = readsex(catname)
    print('Writing to', regname)
    names = s.dtype.names
    if 'X_WORLD' in names and 'Y_WORLD' in names and \
       'A_WORLD' in names and 'B_WORLD' in names and \
       'THETA_WORLD' in names:
        print('Using WORLD coordinates')
        sex_to_DS9reg(regname, s, use_WORLD=1)
    else:
        sex_to_DS9reg(regname, s)


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
