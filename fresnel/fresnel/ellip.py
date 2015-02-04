#Copyright Angel Yanguas-Gil, 2014
"""
ellip: Methods to simulate and manipulate spectroscopic ellipsometry data
"""

import math as m

from fresnel import fresnel
from fresnel.optutils import degtorad, radtodeg

def ellip(angle, waveln, stack, deg=True):
    """Return the ratio of the reflection Fresnel coefficients
    for a stack for a given wavelength and incident angle.

    The stack is formed by a list of tuples (n, k, thickness). In
    the case of the incident and transmitted media, thickness is
    ignored and (n, k) tuples can be passed instead.

    """
    
    rp, rs, tp, ts = fresnel(angle, waveln, stack, deg)
    return rp/rs


def anglestorho(psi, delta):
    """Calculate the ratio of the reflection Fresnel coefficients 
    from the ellipsometric angles psi and delta (in degrees)

    """

    psi = degtorad(psi)
    delta = degtorad(delta)
    tp = m.tan(psi)
    return complex(tp*m.cos(delta), tp*m.sin(delta))


def rhotoangles(rho):
    """Calculate the ellipsometric angles psi and delta from the ratio
    of the fresnel reflection coefficients rp/rs for p and s polarizations

    """

    r = abs(rho)
    psi = m.atan(r)
    psi = radtodeg(psi)
    delta = m.acos(rho.real/r)
    delta = radtodeg(delta)
    if rho.imag/r < 0:
        delta = 2.*m.pi - delta
    return psi, delta


def rhotoepsilon(inc, rho, deg=False):
    """Calculate the pseudodielectric constant given the angle of
    incidence (in radians unless deg flag is set to True)
    and rho, defined as the ratio of the reflection fresnel
    coefficients for p and s polarizations rp/rs

    """

    if deg:
        inc = degtorad(inc)

    try:
        ep = m.tan(inc)**2*(1.-4.*rho*m.sin(inc)**2/(1.+rho)**2)
    except ZeroDivisionError:
        ep = None
    return ep
