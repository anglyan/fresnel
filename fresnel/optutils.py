# Copyright Angel Yanguas-Gil, 2014
"""
optutils: conversion functions

"""

import math as m
import cmath as cm

def radtodeg(a):
    """Radian to degree conversion

    """
    return 180*a/m.pi

def degtorad(a):
    """Degree to radian conversion

    """
    return m.pi*a/180.

def makeindex(n, k):
    """Return the complex refraction index from the real index and
    the attenuation constant
    
    """
    return n + 1J*k

def epton(ep):
    """dielectric constant to (complex) refraction index conversion
    
    """ 
    e1 = ep.real
    e2 = ep.imag
    t2 = m.sqrt(e1*e1+e2*e2)
    n1 = m.sqrt(0.5*(e1+t2))
    n2 = m.sqrt(0.5*(t2-e1))
    return n1+1J*n2

def ntoep(n):
    """(complex) refraction index to dielectric constant conversion
    
    """
    return n*n

def evtonm(en):
    """eV to nm conversion
    
    """
    return 1240./en

def nmtoev(nm):
    """nm to eV conversion
    
    """
    return 1240./nm

def hztoev(hz):
    """Frequency in Hz to energy in eV conversion
    
    """
    return hz*4.136e-15

def evtohz(en):
    """Energy in eV to frequency in Hz conversion
    
    """
    return en/4.136e-15

def cmtoev(k):
    """Wavenumber in cm-1 to energy in eV conversion
    
    """
    return 1.240e-4*k

def evtocm(en):
    """Eenergy in eV to wavenumber in cm-1 conversion
    
    """
    return 8065*en

