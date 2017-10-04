from __future__ import division
from past.utils import old_div
from bpz_tools import *


def function(z, m, nt):
    nz = len(z)
    p_i = ones((nz, nt)) * 1.
    fq = 1.
    fs = 1.
    p_i[1:, 1:nt] = 0.
    ns = sum(p_i[0, 1:nt])
    nq = sum(p_i[:, 0])
    #Normalize relative fractions
    p_i[:, 0] *= (old_div(fq, nq))
    p_i[0, 1:nt] *= (old_div(fs, ns))
    norm = add.reduce(p_i[:nz, :], 0)
    return old_div(p_i[:nz, :], norm[:])
