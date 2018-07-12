'''Function for binning and FB asymmetry calculation'''

#!//////////////////////////////////////////////////////////#
#!//e+e- -> Z/gamma -> mu+mu-
#!//Andreas Papaefstathiou
#!//July 2014
#!//////////////////////////////////////////////////////////#
#! /usr/bin/env python

#python stuff that may or may not be useful
#include "Python.h"
from __future__ import with_statement
from __future__ import division
import vegas
import math, cmath, string, fileinput, pprint
import  os
import sys
from optparse import OptionParser
import random
import cython

#comment out the following if not using matplotlib and numpy
import matplotlib
import mpmath as mp
import numpy as np
import pylab as pl
import scipy
from scipy import interpolate, signal
import matplotlib.font_manager as fm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import sys
from scipy.special import spence as Li
import time
import errno

## some parameters: equivalent to those in MadGraph
pb_convert = 3.894E8 # conversion factor GeV^-2 -> pb
MZ = 91.188    # Z boson mass
GAMMAZ = 2.4414 # Z boson width
alpha = 1/132.507 # alpha QED
Gfermi = 1.16639E-5 # fermi constant
sw2 = 0.222246 # sin^2(weinberg angle)

# e+e- COM energy in GeV
ECM = 90
hats = ECM**2
print "e+e- com energy:", ECM, "GeV"

## define a function that gives the differential cross section

date = time.strftime("%m-%d-%H:%M")
cwd = os.getcwd()
base = os.path.splitext(os.path.basename(__file__))[0]

NEVAL = 50000000
NEV_CON = 9000000
print ('NEVAL = %s' % NEVAL)
print time.strftime("%m-%d-%H:%M")

def dsigma(x):
    # CL and CR
    costh = -1 + (2* x[0])
    phi=  x[1] * 2 * np.pi
    CV = -0.5 + 2 * sw2
    CA = -0.5
    # constants and functions that appear in the differential cross section
    kappa = math.sqrt(2) * Gfermi * MZ**2 / (4 * math.pi * alpha)
    chi1 = kappa * hats * ( hats - MZ**2 ) / (  (hats-MZ**2)**2 + GAMMAZ**2*MZ**2 )
    chi2 = kappa**2 * hats**2 / (  (hats-MZ**2)**2 + GAMMAZ**2*MZ**2 )
    A0 = 1 + 2 * CV**2 * chi1 + (CA**2 + CV**2)**2 * chi2  # see notes
    A1 = 4 * CA**2 * chi1 + 8 * CA**2 * CV**2 * chi2 
    JAC = 4*math.pi
    PREFAC = alpha**2 / (4*hats)*pb_convert * JAC# 2 * pi comes from d-phi integral
    
    return  PREFAC * ( A0 * ( 1 + costh**2 ) + A1 * costh )

integ = vegas.Integrator([[0, 1], [0,1]])
result = integ(dsigma, nitn=10, neval=NEVAL)
print(result.summary())
print('total CS(pb) = %s    Q = %.2f' % (result, result.Q))

summed_res = 0.0
PScosth = []
PSpT=[]
wghts=[]
PSrapidity=[]
Neve = sum(1 for i in integ.random())
jj = 0
for x, wgt in integ.random():
    sys.stdout.write("progress: %d%%   \r" % (float(jj)*100./(Neve)) )
    sys.stdout.flush()
    summed_res += wgt * dsigma(x)
    phi = x[1] * 2 * math.pi
    costh_ii = -1 + (2 * x[0])
    sinphi = math.sin(phi)
    cosphi = math.cos(phi)
    sinth = math.sqrt( 1 - costh_ii**2 )
    pem = [ 0.5 * ECM, 0., 0., 0.5 * ECM ]
    pep = [ 0.5 * ECM, 0., 0., - 0.5 * ECM ]
    pmm = [ 0.5 * ECM, 0.5 * ECM * sinth * cosphi, 0.5 * ECM * sinth * sinphi, 0.5 * ECM * costh_ii ]
    pmp = [ 0.5 * ECM, - 0.5 * ECM * sinth * cosphi, - 0.5 * ECM * sinth * sinphi, - 0.5 * ECM * costh_ii ]
    pTmu = (pmm[1]**2 + pmm[2]**2)**(0.5)
    calc_rapidity = 0.5 * np.log((pmm[0] + pmm[3])/(pmm[0] - pmm[3]))

    PScosth.append(costh_ii)
    PSpT.append(pTmu)
    wghts.append(dsigma(x)*wgt)
    PSrapidity.append(calc_rapidity)
    jj +=1

print ('summed sigma(pb) = %s' % summed_res)

'----=========================================================================----'
'----==========================      ANALYTICAL     ==========================----'
'----=========================================================================----'
def Analytic(x1, x0):
# CL and CR
    CV = -0.5 + 2 * sw2
    CA = -0.5
    # constants and functions that appear in the differential cross section
    kappa = math.sqrt(2) * Gfermi * MZ**2 / (4 * math.pi * alpha)
    chi1 = kappa * hats * ( hats - MZ**2 ) / (  (hats-MZ**2)**2 + GAMMAZ**2*MZ**2 )
    chi2 = kappa**2 * hats**2 / (  (hats-MZ**2)**2 + GAMMAZ**2*MZ**2 )
    A0 = 1 + 2 * CV**2 * chi1 + (CA**2 + CV**2)**2 * chi2  # see notes
    A1 = 4 * CA**2 * chi1 + 8 * CA**2 * CV**2 * chi2
    ANALYTIC = math.pi*alpha**2*(A0*(x1-x0)+A0*(x1**3-x0**3)/3+A1*(x1**2-x0**2)/2)/2/hats * pb_convert
    return ANALYTIC

print 'total analytic CS(pb) = %s' % Analytic(1, -1)

'----=========================================================================----'
'----==========================      SAVE FILES     ==========================----'
'----=========================================================================----'
 #save array
if NEVAL > NEV_CON:    
    ee_array = [PScosth, PSpT, PSrapidity, wghts]
     #make directory path

    Save_path = cwd
    SaveStr = Save_path + '/ee_data'         
    np.save(SaveStr, ee_array)

     #save text
    with open(os.path.join(Save_path, 'details.txt'), 'a') as myfile:
        myfile.write('\n\n || total CS(pb) = %s || total analytic CS(pb) = %s || NEVAL = %s ||\n ' % (summed_res, Analytic(1,-1), NEVAL))

print time.strftime("%m-%d-%H:%M")



            
