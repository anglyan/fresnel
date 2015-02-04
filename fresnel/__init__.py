#Copyright Angel Yanguas-Gil, 2014

"""Python Fresnel package

Fresnel is a python package that calculates coherent transmission
and reflexion of light from a multilayer structure.

Modules_

    base: contains the methods required to determine the reflection and
    transmission coefficients of interfaces and stacks
    optutils: conversion functions
    ellip: methods to model the spectroscopic ellipsometry response of
    a multilayer system

Methods_

    fresnel
    evtonm
    nmtoev
    ntoep
    epton

"""

from fresnel.base import fresnel, fresnel2, rs, rp
from fresnel.optutils import evtonm, nmtoev, ntoep, epton

