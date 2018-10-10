import string
import numpy

""" A collection of scripts to read and write very simple ascii tables.
    Collected by F.Menanteau.
"""

#Read/write headers


def get_header(file):
    """ Returns a string containing all the lines
    at the top of a file which start by '#'"""
    buffer = ''
    for line in open(file).readlines():
        if line[0] == '#':
            buffer = buffer + line
        else:
            break
    return buffer


def put_header(file, text):
    """Adds text (starting by '#' and ending by '\n')
    to the top of a file."""
    if len(text) == 0:
        return
    if text[0] != '#':
        text = '#' + text
    if text[-1] != '\n':
        text = text + '\n'
    buffer = text + open(file).read()
    open(file, 'w').write(buffer)

#Files containing strings


def get_str(file, cols=0, nrows='all', sep=None):
    """
        Reads strings from a file
        Usage:
         x,y,z=get_str('myfile.cat',(0,1,2))
        x,y,z are returned as string lists
    """
    # Modified to be feed a buffer as well as a file
    # F. Menanteau

    if isinstance(cols, int):
        cols = (cols, )
        nvar = 1
    else:
        nvar = len(cols)
    lista = []
    for i in range(nvar):
        lista.append([])

    if type(file) is list:
        buffer = file
        #print "# Passing a buffer"
    else:
        buffer = open(file).readlines()

    if nrows == 'all':
        nrows = len(buffer)
    counter = 0
    for lines in buffer:
        if counter >= nrows:
            break
        if sep:
            pieces = lines.split(sep)
        else:
            pieces = lines.split()

        if len(pieces) == 0:
            continue
        if pieces[0][0] == '#':
            continue
        for j in range(nvar):
            lista[j].append(pieces[cols[j]])
        counter = counter + 1
    if nvar == 1:
        return lista[0]
    else:
        return tuple(lista)


def put_str(file, tupla):
    """ Writes tuple of string lists to a file
        Usage:
      put_str(file,(x,y,z))
    """
    if not isinstance(tupla, tuple):
        raise 'Need a tuple of variables'

    f = open(file, 'w')

    for i in range(1, len(tupla)):
        if len(tupla[i]) != len(tupla[0]):
            raise 'Variable lists have different length'
    for i in range(len(tupla[0])):
        cosas = []
        for j in range(len(tupla)):
            cosas.append(str(tupla[j][i]))
        f.write('\n'.join(cosas))
    f.close()

#Files containing data


def get_data(file, cols=0, nrows='all', sep=None):
    """ Returns data in the columns defined by the tuple
    (or single integer) cols as a tuple of float arrays
    (or a single float array)"""

    if isinstance(cols, int):
        cols = (cols, )
        nvar = 1
    else:
        nvar = len(cols)

    data = get_str(file, cols, nrows, sep=sep)

    if nvar == 1:
        return numpy.array(list(map(float, data)))
    else:
        data = list(data)
        for j in range(nvar):
            data[j] = numpy.array(list(map(float, data[j])))
        return tuple(data)


def put_data(file, variables, header='', format='', append='no'):
    """ Writes tuple of float variables to a file
        Usage:
      put_data(file,(x,y,z),header,format)
    where header is any string
        and format is a string of the type:
           '%f %f %i '
    """
    if not isinstance(variables, tuple):
        raise 'Need a tuple of variables'
    if format == '': format = '%.8e  ' * len(variables)
    if append == 'yes': f = open(file, 'a')
    else: f = open(file, 'w')
    if header != "":
        if header[0] != '#': header = '#' + header
        if header[-1] != '\n': header = header + '\n'
        f.write(header)
    for i in range(len(variables[0])):
        cosas = []
        for j in range(len(variables)):
            cosas.append(variables[j][i])
        line = format % tuple(cosas)
        f.write(line + '\n')
    f.close()

#Files containing strings

# F. Menanteau, reads all or some columns


def rcols(file, cols=None, nrows='all'):
    """ Returns data in the columns defined by the tuple
    (or single integer) cols as a tuple of float arrays
    (or a single float array)"""

    if cols is None:
        nvar = 0
    elif isinstance(cols, int):
        cols = (cols, )
        nvar = 1
    else:
        nvar = len(cols)

    data = get_string(file, cols, nrows)

    if nvar == 1:
        return numpy.array(list(map(float, data)))
    else:
        data = list(data)
        nvar = len(data)
        for j in range(nvar):
            data[j] = numpy.array(list(map(float, data[j])))
        return tuple(data)


def get_string(file, cols=None, nrows='all', buffer=None):
    """
        Reads strings from a file
        Usage:
         x,y,z=get_str('myfile.cat',(0,1,2))
        x,y,z are returned as string lists

        Modified to read from buffer, F. Menanteau

    """
    nvar = None

    if (buffer):
        buffer = file
    else:
        buffer = open(file).readlines()
    if nrows == 'all':
        nrows = len(buffer)
    counter = 0
    for lines in buffer:

        if counter >= nrows:
            break
        pieces = string.split(lines)
        if len(pieces) == 0:
            continue
        if pieces[0][0] == '#':
            continue

    # Decide how many columns to read
        if nvar is None:

            if cols is None:
                nvar = len(pieces)
                cols = tuple(range(nvar))
            elif isinstance(cols, int):
                cols = (cols, )
                nvar = 1
            else:
                nvar = len(cols)
            lista = []
            for i in range(nvar):
                lista.append([])

        for j in range(nvar):
            lista[j].append(pieces[cols[j]])
        counter = counter + 1

    if nvar == 1:
        return lista[0]
    else:
        return tuple(lista)


def get_datarray(file, cols=0, nrows='all', buffer=None):
    """ Returns data in the columns defined by the tuple
    (or single integer) cols as a tuple of float arrays
    (or a single float array)

    Modified to read from buffer, F. Menanteau"""

    if isinstance(cols, int):
        cols = (cols, )
        nvar = 1
    else:
        nvar = len(cols)

    data = get_string(file, cols, nrows, buffer)

    if nvar == 1:
        return numpy.array(list(map(float, data)))
    else:
        data = list(data)
        for j in range(nvar):
            data[j] = numpy.array(list(map(float, data[j])))
        return tuple(data)

#Read/write 2D arrays
# Added from useful.py


def get_2Darray(file, cols='all', nrows='all', verbose='no'):
    """Read the data on the defined columns of a file
    to an 2 array
    Usage:
    x=get_2Darray(file)
    x=get_2Darray(file,range(len(p))
    x=get_2Darray(file,range(0,10,2),nrows=5000)
    Returns x(nrows,ncols)
    """
    if cols == 'all':
        #Get the number of columns in the file
        for line in open(file).readlines():
            pieces = string.split(line)
            if len(pieces) == 0: continue
            if line[0] == '#': continue
            nc = len(pieces)
            cols = list(range(nc))
            if verbose == 'yes': print('cols=', cols)
            break
    else:
        nc = len(cols)

    lista = get_data(file, cols, nrows)
    nl = len(lista[0])
    x = numpy.zeros((nl, nc), type='Float64')  # Float64 to avoid warning msgs.
    for i in range(nc):
        x[:, i] = lista[i]

    return x


def put_2Darray(file, array, header='', format='', append='no'):
    """ Writes a 2D array to a file, where the first
    index changes along the lines and the second along
    the columns
    Usage: put_2Darray(file,a,header,format)
    where header is any string
        and format is a string of the type:
           '%f %f %i '
    """
    lista = []
    for i in range(array.shape[1]):
        lista.append(array[:, i])
    lista = tuple(lista)
    put_data(file, lista, header, format, append)
