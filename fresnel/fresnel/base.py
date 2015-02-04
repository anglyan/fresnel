#Copyright Angel Yanguas-Gil, 2014
"""
base: Fresnel conefficients and reflection and transmission
of polarized light from interfaces and multilayers determined using
the transfer matrix method.

"""
from __future__ import division

import cmath as cm
import math as m

try:
    from optutils import makeindex, radtodeg, degtorad
except ImportError:
    def makeindex(n0, n1):
        return n0 + 1J*n1

    def radtodeg(a):
        return m.pi*a/180.

    def degtorad(a):
        return a*180./m.pi


def fresnel(t0, wl, layers, deg=False):
    """
    Determine the fresnel coefficients of a multilayer.

    Arguments:
    t0: incident angle (in radians or degrees if deg=True)
    wl: wavelength (same units as thickness)
    layers: list [(n, k), (n, k, thickness),..., (n, k)]
    deg: flag indicating whether the incidence angle is in degrees.
    
    Return the tuple (rp, rs, tp, ts) with the reflection
    and transmission Fresnel coefficients for the p and s
    polarizations.

    """

    if deg == True:
        t0 = t0*m.pi/180.

#    if wl <= 0 or t0 >= 0.5*m.pi or t0 <= -0.5*m.pi:
#        raise ValueError("Parameter out of range")

    if len(layers) == 1:
        return (0.+0J, 0.+0J, 1.+0J, 1.+0J)
    elif len(layers) == 2:
        return fresnel2(t0, makeindex(layers[0][0], layers[0][1]),
            makeindex(layers[1][0], layers[1][1]))
    else:
        return _fresneln(t0, wl, layers)


def rp(t0, N0, N1, deg=False):
    """
    Calculate the reflection coefficient for p polarization
    
    Arguments:
    t0: incident angle (in radians or degrees if deg=True)
    N0: complex refraction index of the incident medium (n+ik)
    N1: complex refraction index of the second medium (n+ik)
    deg: flag controlling the incidence angle units


    Return rp, rs, tp, ts for the p and s polarizations.
    
    """

    if deg == True:
        t0 = t0*m.pi/180.
    ct0 = cm.cos(t0)
    st0 = cm.sin(t0)
    st1 = N0/N1*st0
    ct1 = cm.sqrt(1-st1*st1)
    rp = (N1*ct0 - N0*ct1)/(N1*ct0 + N0*ct1)
    return rp
    
def rs(t0, N0, N1, deg=False):
    """
    Calculate the reflection coefficient for p polarization
    
    Arguments:
    t0: incident angle (in radians or degrees if deg=True)
    N0: complex refraction index of the incident medium (n+ik)
    N1: complex refraction index of the second medium (n+ik)
    deg: flag controlling the incidence angle units


    Return rp, rs, tp, ts for the p and s polarizations.
    
    """

    if deg == True:
        t0 = t0*m.pi/180.
    ct0 = cm.cos(t0)
    st0 = cm.sin(t0)
    st1 = N0/N1*st0
    ct1 = cm.sqrt(1-st1*st1)
    rs = (N0*ct0 - N1*ct1)/(N0*ct0 + N1*ct1)
    return rs

def fresnel2(t0, N0, N1, deg=False):
    """
    Calculate the Fresnel coefficients for an interface.
    
    Arguments:
    t0: incident angle (in radians or degrees if deg=True)
    N0: complex refraction index of the incident medium (n+ik)
    N1: complex refraction index of the second medium (n+ik)
    deg: flag controlling the incidence angle units


    Return rp, rs, tp, ts for the p and s polarizations.
    
    """

    if deg == True:
        t0 = t0*m.pi/180.
    ct0 = cm.cos(t0)
    st0 = cm.sin(t0)
    st1 = N0/N1*st0
    ct1 = cm.sqrt(1-st1*st1)
    rp = (N1*ct0 - N0*ct1)/(N1*ct0 + N0*ct1)
    rs = (N0*ct0 - N1*ct1)/(N0*ct0 + N1*ct1)
    tp = (2.*N0*ct0)/(N1*ct0 + N0*ct1)
    ts = (2.*N0*ct0)/(N0*ct0 + N1*ct1)
    return rp, rs, tp, ts


def _fresneln(t0, wl, layers):
    """
    
    Determine the transmission and reflection coefficients of 
    fresnel coefficients of a multilayer for two or
    more interfaces.

    Arguments:
    t0: incident angle (in radians or degrees if deg=True)
    wl: wavelength (same units as thickness)
    layers: list [(n, k), (n, k, thickness),..., (n, k)]
    
    It returns the tuple (rp, rs, tp, ts) with the reflection
    and transmission Fresnel coefficients for the p and s
    polarizations.
    
    """

    N0 = makeindex(layers[0][0], layers[0][1])
    NS = makeindex(layers[-1][0], layers[-1][1])
    layers = layers[1:-1]
    
    mlp = []
    mls = []
    Nold = N0
    told = t0
    for layer in layers:
        Nnew = makeindex(layer[0], layer[1])
        mp, ms = _matinterface(Nold, Nnew, told)
        tnew = cm.asin(Nold*cm.sin(told)/Nnew)
        dnorm = layer[2]/wl*cm.cos(tnew)
        ml = _matlayer(Nnew, dnorm)
        mlp.append(mp)
        mls.append(ms)
        mlp.append(ml)
        mls.append(ml)
        Nold = Nnew
        told = tnew

    mp, ms = _matinterface(Nold, NS, told)
    mlp.append(mp)
    mls.append(ms)

    mp = _nmatmult(mlp)
    ms = _nmatmult(mls)

    rp = mp[2]/mp[0]
    rs = ms[2]/ms[0]
    tp = 1./mp[0]
    ts = 1./ms[0]

    return rp, rs, tp, ts


def fresnelreverse(rp, rs, tp, ts):
    """
    Calculate the inverse fresnel coefficients for an interface
    given the fresnel coefficients in one direction

    Arguments:
    rp, rs, tp, ts: reflection and transmission fresnel coefficients
    for the p and s polarizations.

    Return the inverse Fresnel coefficients rp2, rs2, tp2, ts2
    
    """
    
    rp2 =- rp
    rs2 =- rs
    tp2 = (1.-rp*rp)/tp
    ts2 = (1.-rs*rs)/ts
    return rp2, rs2, tp2, ts2


def _matinterface(N0, Nl, t0):
    """Return the matrix for an interface between media with complex 
    refractive indices N0 and N1 and incidence angle t0.
    
    """
    
    rp, rs, tp, ts = fresnel2(t0, N0, Nl)
    m1p = 1./tp
    m1s = 1./ts
    m2p = rp/tp
    m2s = rs/ts
    mp = (m1p, m2p, m2p, m1p)
    ms = (m1s, m2s, m2s, m1s)
    return mp, ms

def _matlayer(N, dn):
    """Return the matrix for a layer of refractive index N and 
    thickness (normalized to the vacuum wavelength) dn.
    
    """

    beta = 2.*m.pi*dn*N
    b1 = cm.exp(-1J*beta)
    b2 = cm.exp(1J*beta)
    return (b1, 0., 0., b2)

def _nmatmult(ml):
    """Multiply a list of 2x2 matrices coded as 4-tuples."""
    return reduce(_matmult, ml)

def _matmult(m1, m2):
    """Multiply two 2x2 matrices presented as 4-tuples."""
    m11 = m1[0]*m2[0]+m1[1]*m2[2]
    m12 = m1[0]*m2[1]+m1[1]*m2[3]
    m21 = m1[2]*m2[0]+m1[3]*m2[2]
    m22 = m1[2]*m2[1]+m1[3]*m2[3]
    return m11, m12, m21, m22

if __name__ == '__main__':
    mats = [(1., 0.), (2., 3., 100.), (1., 0.)]
    rp, rs, tp, ts = fresnel(0., 400., mats)
    print tp
    print ts

