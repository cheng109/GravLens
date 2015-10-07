#!/usr/bin/env python

import sys
from math import *

try:
  H0 = 72                         # Hubble constant
  WM = 0.26                        # Omega(matter)
  WV =  0 #1.0 - WM - 0.4165/(H0*H0)  # Omega(vacuum) or lambda
  z = 0.445

  WR = 0.        # Omega(radiation)
  WK = 0.        # Omega curvaturve = 1-Omega(total)
  c = 299792.458 # velocity of light in km/sec


  DCMR = 0.0     # comoving radial distance in units of c/H0
  DCMR_Mpc = 0.0
  DCMR_Gyr = 0.0
  DA = 0.0       # angular size distance
  DA_Mpc = 0.0
  DA_Gyr = 0.0
  kpc_DA = 0.0


  h = H0/100.
  WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
  WK = 1-WM-WR-WV
  az = 1.0/(1+1.0*z)

  n=1000         # number of points in integrals

# do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
  for i in range(n):
    a = az+(1-az)*(i+0.5)/n
    adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))

    DCMR = DCMR + 1./(a*adot)


  DCMR = (1.-az)*DCMR/n

  ratio = 1.00
  x = sqrt(abs(WK))*DCMR
  if x > 0.1:
    if WK > 0:
      ratio =  0.5*(exp(x)-exp(-x))/x
    else:
      ratio = sin(x)/x
  else:
    y = x*x
    if WK < 0: y = -y
    ratio = 1. + y/6. + y*y/120.
  DCMT = ratio*DCMR
  DA = az*DCMT
  DA_Mpc = (c/H0)*DA


except:
    print "hello"



print DA_Mpc