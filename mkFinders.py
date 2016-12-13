import numpy as np
import pyds9
from astLib import astCoords
from datetime import date

# load the data
data = np.genfromtxt(
    'PSZ2_unconfirmed_catalog_4NOAO_2016A_newSwift - PSZ2_unconfirmed_catalog_4NOAO_.csv',
    delimiter=',',
    names=True,
    dtype=None)

# set up the ds9 enviroment
params = [
    'view colorbar no', 'cmap invert yes', 'grid yes', 'grid view grid no',
    'grid view title no', 'grid type publication', 'grid type axes exterior',
    'grid type numerics exterior', 'view info no', 'view panner no',
    'view magnifier no', 'view buttons no', 'page setup orientation portrait',
    'page setup pagesize letter'
]

# start the ds9 session
d = pyds9.DS9()

for param in params:
    d.set(param)

d.set('dssstsci survey poss2ukstu_red')
d.set('dssstsci size 45.0 45.0 arcmin')

for i, _ in enumerate(data['RA']):
    d.set('dssstsci coord %f %f' % (data['RA'][i], data['Dec'][i]))
    d.set('orient x')
    d.set('zoom to fit')
    d.set('regions', "image; text 1450 2700 # text={%s Big Field - %s} \
        font={times 18 bold} color = black" % (data['Name'][i], date.today()))
    d.set('regions',
          'compass(2560,2350,30) #compass=fk5 {N} {E} 1 1 color=black')

    # draw a circle at the planck position
    cir = "fk5;circle(%f,%f,3') #color=black" % (data['RA'][i], data['Dec'][i])
    d.set('regions', cir)

    # Upper left
    newRA, newDec = astCoords.shiftRADec(data['RA'][i], data['Dec'][i],
                                         -8 * 60, 8 * 60)
    box = "fk5; box(%f, %f, 28', 28', 0) #color=green" % (newRA, newDec)
    d.set('regions', box)

    # Upper right
    newRA, newDec = astCoords.shiftRADec(data['RA'][i], data['Dec'][i], 7 * 60,
                                         7 * 60)
    box = "fk5; box(%f, %f, 28', 28', 0) #color=blue" % (newRA, newDec)
    d.set('regions', box)

    # lower right
    newRA, newDec = astCoords.shiftRADec(data['RA'][i], data['Dec'][i],
                                         -7 * 60, -7 * 60)
    box = "fk5; box(%f, %f, 28', 28', 0) #color=black" % (newRA, newDec)
    d.set('regions', box)

    # lower left
    newRA, newDec = astCoords.shiftRADec(data['RA'][i], data['Dec'][i], 8 * 60,
                                         -8 * 60)
    box = "fk5; box(%f, %f, 28', 28', 0) #color=red" % (newRA, newDec)
    d.set('regions', box)

    # # 4 pointing
    # box = "fk5; box(%f, %f, 14', 14', 0) #color=red" % (data['RA'][i],
    #                                                     data['Dec'][i])
    # d.set('regions', box)
    #
    # # 2 pointing
    # newRA, newDec = astCoords.shiftRADec(data['RA'][i], data['Dec'][i], -14*60, 0)
    # box = "fk5; box(%f, %f, 14', 14', 0) #color=green" % (newRA,newDec)
    # d.set('regions', box)
    #
    # newRA, newDec = astCoords.shiftRADec(data['RA'][i], data['Dec'][i], 14*60, 0)
    # box = "fk5; box(%f, %f, 14', 14', 0) #color=green" % (newRA,newDec)
    # d.set('regions', box)
    #
    # newRA, newDec = astCoords.shiftRADec(data['RA'][i], data['Dec'][i], 0, -14*60)
    # box = "fk5; box(%f, %f, 14', 14', 0) #color=green" % (newRA,newDec)
    # d.set('regions', box)
    #
    # newRA, newDec = astCoords.shiftRADec(data['RA'][i], data['Dec'][i], 0, 14*60)
    # box = "fk5; box(%f, %f, 14', 14', 0) #color=green" % (newRA,newDec)
    # d.set('regions', box)
    #
    # # single pointing
    # newRA, newDec = astCoords.shiftRADec(data['RA'][i], data['Dec'][i], -14*60,
    #                                      -14*60)
    # box = "fk5; box(%f, %f, 14', 14', 0) #color=magenta" % (newRA,newDec)
    # d.set('regions', box)
    #
    # newRA, newDec = astCoords.shiftRADec(data['RA'][i], data['Dec'][i], -14*60,
    #                                      14*60)
    # box = "fk5; box(%f, %f, 14', 14', 0) #color=magenta" % (newRA,newDec)
    # d.set('regions', box)
    #
    # newRA, newDec = astCoords.shiftRADec(data['RA'][i], data['Dec'][i], 14*60, -14*60)
    # box = "fk5; box(%f, %f, 14', 14', 0) #color=magenta" % (newRA,newDec)
    # d.set('regions', box)
    #
    # newRA, newDec = astCoords.shiftRADec(data['RA'][i], data['Dec'][i], 14*60, 14*60)
    # box = "fk5; box(%f, %f, 14', 14', 0) #color=magenta" % (newRA,newDec)
    # d.set('regions', box)

    d.set('raise')
    d.set('regions select none')
    d.set('print resolution 75')
    d.set('print level 2')
    d.set('print palette gray')
    d.set('print destination file')
    d.set('print filename ./images/%s_big.ps' % data['Name'][i])
    d.set('print')

    # save the regions
    d.set('regions system wcs')
    d.set('regions sky fk5')
    d.set('regions save ./regions/%s_big.reg' % data['Name'][i])
    d.set('frame clear')
