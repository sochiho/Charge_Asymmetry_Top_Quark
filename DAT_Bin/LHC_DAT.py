'''Function for obtaining TEVATRON data, collecting the x variables plugged into VEGAS'''
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
import cPickle as pickle
import shutil
import resource
import gc

sys.path.append('/home/chiho/lhapdf/lib/python2.7/site-packages')

import lhapdf

print '\n'
print '----====================================================----'
print 'Exercise 2 (Advanced Scientific Computing Workshop, July 2014)'
print 'pp_bar --> g --> t t_bar'
print '----====================================================----'
print '\n'
    
## some parameters
pb_convert = 3.894E8 # conversion factor GeV^-2 -> pb
MZ = 91.188    # Z boson mass
GAMMAZ = 2.4414 # Z boson width
alpha = 1/132.507 # alpha QED
Gfermi = 1.16639E-5 # fermi constant
sw2 = 0.222246 # sin^2(weinberg angle)

#Define variables in Kuhn paper
M_top = 173.1 # Mass of top quark GeV
E_cut = 0.01 * M_top

# we also need the "ranges" (i.e. x2 - x1)
deltath = 2

def mem():
    print('Memory usage         : % 8.3f GB' % (
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/10**6))

'----=========================================================================----'
'----====================     PDF SYMMETRIC FUNCTION    ======================----'
'----=========================================================================----'

def fqfqbar_symmetric(x1, x2):
    '''Function for the qqbar symmetric PDFs'''
    mu = M_top
    ff = 0
    quark = [1, 2, 3, 4]
    for ii in range(len(quark)):
        i = quark[ii]
        ff += (p.xfxQ(i, x1, mu) *  p.xfxQ(-i, x2, mu)) + (p.xfxQ(-i, x1, mu) *  p.xfxQ(i, x2, mu))
    return ff /(x1 * x2)

def fgfg_symmetric(x1, x2):
    '''Function for the gg symmetric PDFs'''
    mu = M_top
    ff = 0
    ff += (p.xfxQ(21, x1, mu) *  p.xfxQ(21, x2, mu)) # + (p.xfxQ(-i, x1, mu) *  p.xfxQ(-i, x2, mu))
    return ff /(x1 * x2)

'----=========================================================================----'
'----===========================    LO FUNCTIONS    ==========================----'
'----=========================================================================----'
##Calulate Leading Order terms
## define a function that gives the differential cross section at LO
## as a function of the parameters
def Top_LO_qq(x):
    #x1 and x2 values
    x1 = tau0 + (1 -tau0) * x[0]
    x2 = tau0/x1 +(1-tau0/x1)*x[1]
    hats = x1*x2*s
    if hats < 4*M_top**2:
        print 'hats < 4M_top'
    xm = M_top **2 / hats
    Q = math.sqrt(hats)
    
    #random costh value
    costh = -1 + deltath*x[2]
    sinth = math.sqrt(1- costh**2)
    beta = math.sqrt(1 - (4 * xm))
    c = beta * costh
    
    E3 = 0.5
    k3 = math.sqrt(E3**2 -xm)
    y3 = np.log((E3 + k3*costh)/(E3 - k3*costh))/2 #tt frame
    y4 = -y3    
    #in Lab frame 
    y3_lab = y3 + np.log(x1/x2)/2
    y4_lab = -y3 + np.log(x1/x2)/2

    #psuedo-rapidity
    y3p = np.log((1 + costh)/(1 - costh))/2
    y3p_lab = y3p + np.log(x1/x2)/2
    y4p_lab = -y3p + np.log(x1/x2)/2

    #ABS FRAME is symmetric anyway so gives same answer as just integrating
    if FRAME == 'ABS_y':
        CON = abs(y3_lab) - abs(y4_lab) #abs(y3) - abs(-y3)

    elif FRAME == 'PSE_RAP':
        CON = abs(y3p_lab) - abs(y4p_lab)
        
    elif FRAME == 'OUT_CUT':
        CON1 = abs(y3_lab) - Y_CUT
        CON2 = abs(y4_lab) - Y_CUT
        
    elif FRAME == 'IN_CUT':
        CON1 = Y_CUT - abs(y4_lab)
        CON2 = Y_CUT -abs(y3_lab)
        
    #Impose further cuts on y3
    elif FRAME == 'OUT_CUTC':
        CON1 = 0
        CON2 = 0
        if abs(y3_lab) < 2.5:
            CON1 = abs(y3_lab) - Y_CUT
        if abs(y4_lab) < 2.5: 
            CON2 = abs(y4_lab) - Y_CUT
            
    elif FRAME == 'Y_CUT':
        CON1 = 0
        CON2 = 0
        if (y3_lab + y4_lab)/2 > Y_CUT:
            CON1 = y3_lab - y4_lab
            CON2 = y4_lab - y3_lab

    elif FRAME == 'Y_CUTP':
        pT = k3 * sinth * Q * 0.5
        CON1 = 0
        CON2 = 0
        if (y3_lab + y4_lab)/2 > Y_CUT and pT < 20:
            CON1 = y3_lab - y4_lab
            CON2 = y4_lab - y3_lab

    elif FRAME == 'Y_CUTQ':
        CON1 = 0
        CON2 = 0
        if (y3_lab + y4_lab)/2 > Y_CUT and Q > 450:
            CON1 = y3_lab - y4_lab
            CON2 = y4_lab - y3_lab
            
    rf = 0
    rb = 0
                
    if FRAME in FRAME_lB or FRAME in FRAME_lC:
        CON = CON1

    if CON > 0:
        rf += (1 + c**2 + (4 * xm))

    if FRAME in FRAME_lB or FRAME in FRAME_lC:
        CON = -CON2    
    if CON < 0:
        rb += (1 + c**2 + (4 * xm))

    w_ii = rf + rb
    JAC = (1 -tau0) * (1-tau0/x1) * deltath
    PREFAC = alpha_QCD**2 * math.pi * beta * pb_convert /(9 * hats) 
    return  PREFAC * JAC * w_ii * fqfqbar_symmetric(x1, x2)

def Top_LO_gg(x):
    x1 = tau0 + (1 -tau0) * x[0]
    x2 = tau0/x1 +(1-tau0/x1)*x[1]
    hats = x1*x2*s
    xm = M_top **2 / hats
    if hats < 4*M_top**2:
        print 'hats < 4M_top'
    Q = math.sqrt(hats)
        
    #random costh value
    costh = -1 + deltath*x[2]
    sinth = math.sqrt(1- costh**2)
    beta = math.sqrt(1 - (4 * xm))
    c = beta * costh

    E3 = 0.5
    k3 = math.sqrt(E3**2 -xm)
    y3 = np.log((E3 + k3*costh)/(E3 - k3*costh))/2 #tt frame
    y4 = -y3
    
    #in Lab frame 
    y3_lab = y3 + np.log(x1/x2)/2
    y4_lab = -y3 + np.log(x1/x2)/2
    #psuedo-rapidity
    y3p = np.log((1 + costh)/(1 - costh))/2
    y3p_lab = y3p + np.log(x1/x2)/2
    y4p_lab = -y3p + np.log(x1/x2)/2
    
    #ABS FRAME is symmetric anyway so gives same answer as just integrating
    if FRAME == 'ABS_y':
        CON = abs(y3_lab) - abs(y4_lab) #abs(y3) - abs(-y3)

    elif FRAME == 'PSE_RAP':
        CON = abs(y3p_lab) - abs(y4p_lab)

    elif FRAME == 'OUT_CUT':
        CON1 = abs(y3_lab) - Y_CUT
        CON2 = abs(y4_lab) - Y_CUT

    elif FRAME == 'IN_CUT':
        CON1 = Y_CUT - abs(y4_lab)
        CON2 = Y_CUT -abs(y3_lab)

    #Impose further cuts on y3
    elif FRAME == 'OUT_CUTC':
        CON1 = 0
        CON2 = 0
        if abs(y3_lab) < 2.5:
            CON1 = abs(y3_lab) - Y_CUT
        if abs(y4_lab) < 2.5: 
            CON2 = abs(y4_lab) - Y_CUT
            
    elif FRAME == 'Y_CUT':
        CON1 = 0
        CON2 = 0
        if (y3_lab + y4_lab)/2 > Y_CUT:
            CON1 = y3_lab - y4_lab
            CON2 = y4_lab - y3_lab

    elif FRAME == 'Y_CUTP':
        pT = k3 * sinth * Q * 0.5
        CON1 = 0
        CON2 = 0
        if (y3_lab + y4_lab)/2 > Y_CUT and pT < 20:
            CON1 = y3_lab - y4_lab
            CON2 = y4_lab - y3_lab

    elif FRAME == 'Y_CUTQ':
        CON1 = 0
        CON2 = 0
        if (y3_lab + y4_lab)/2 > Y_CUT and Q > 450:
            CON1 = y3_lab - y4_lab
            CON2 = y4_lab - y3_lab 

    rf = 0
    rb = 0
    
    if FRAME in FRAME_lB or FRAME in FRAME_lC:
        CON = CON1
    if CON > 0:
        rf += (1/(3 * (1-c**2)) - 3/16) * (1 + c**2 + 8*xm - (32*xm**2)/(1-c**2))

    if FRAME in FRAME_lB or FRAME in FRAME_lC:
        CON = -CON2    
    if CON < 0:
        rb += (1/(3 * (1-c**2)) - 3/16) * (1 + c**2 + 8*xm - (32*xm**2)/(1-c**2))

    w_ii = rf + rb
    JAC = (1 -tau0) * (1-tau0/x1) * deltath
    PREFAC = alpha_QCD**2 * math.pi * beta * pb_convert /(2 * hats)
    return  PREFAC * JAC * w_ii * fgfg_symmetric(x1, x2)

'----=========================================================================----'
'----====================     PDF ASYMMETRIC FUNCTION   ======================----'
'----=========================================================================----'

def fqfqbar_asymmetric(x1, x2):
    '''Function for the qqbar asymmetric PDFs'''
    mu = M_top
    ff = 0
    quark = [1, 2, 3, 4]
    for ii in range(len(quark)):
        i = quark[ii]
        ff += (p.xfxQ(i, x1, mu) *  p.xfxQ(-i, x2, mu)) - (p.xfxQ(-i, x1, mu) *  p.xfxQ(i, x2, mu))
    return ff /(x1 * x2)

'----=========================================================================----'
'----======================      SOFT + VIRTUAL     ==========================----'
'----=========================================================================----'

def qqV_B(c, hats):
    '''Function B(c) for soft + virtual contribution'''
    xm = M_top**2 / hats
    b = math.sqrt(1 - (4 * xm))
    
    qqV_B = np.log((1-c)/2) * (1-c**2-8*xm)/(1-c-2*xm) + (c+2*xm) * (2*Li(2*xm/(1-c)) - np.log((1-c)/2)**2) \
            + xm*np.log(xm)*4*c*(2-c**2-7*xm)/b**2/((1-2*xm)**2-c**2) + np.log(xm)**2*c/2 \
            - (np.log((1-b)/(1+b))**2 + 4*Li(1+((1-b)/(1+b))) + np.pi**2/3) * (1-8*xm+8*xm**2)*c/2/b**3 \
            -c*np.pi**2/6
    return qqV_B

def qqV_D(c, hats):
    '''Function D(c) for soft + virtual contribution'''
    xm = M_top**2 / hats
    b = math.sqrt(1 - (4 * xm))
    xx = (1-c)/math.sqrt(2*(1-c-2*xm))
    yy = (1-b+math.sqrt(2*(1-c-2*xm)))/2
    
    qq_V_D = 2*(Li(1+xx/(1-yy)+0j).real- Li(1-(1-xx)/(1-yy)+0j).real - Li(1-(1+xx)/yy+0j).real\
                + Li(1-xx/yy+0j).real)+ np.log(abs(yy/(1-yy)))**2 - Li(1-xx**2+0j).real \
                + np.log(xx**2)**2/2 -np.log(xx**2)*np.log(1-xx**2)
    
    return qq_V_D

def Top_Asym_qqV(x):
    '''Function for soft+virtual to be integrated'''
    x1 = tau0 + (1 -tau0) * x[0]
    x2 = tau0/x1 +(1-tau0/x1)*x[1]
    hats = x1*x2*s
    xm = M_top**2/hats
    Q = math.sqrt(hats)
    b = math.sqrt(1-4*xm)
    w = E_cut/math.sqrt(hats)

    costh = -1 + deltath*x[2]
    sinth = math.sqrt(1- costh**2)
    c = b*costh

    E3 = 0.5
    k3 = math.sqrt(E3**2 -xm)
    y3 = np.log((E3 + k3*costh)/(E3 - k3*costh))/2 #tt frame
    y4 = -y3
        
    #in Lab frame 
    y3_lab = y3 + np.log(x1/x2)/2
    y4_lab = -y3 + np.log(x1/x2)/2
    #psuedo-rapidity
    y3p = np.log((1 + costh)/(1 - costh))/2
    y3p_lab = y3p + np.log(x1/x2)/2
    y4p_lab = -y3p + np.log(x1/x2)/2
    
    if FRAME == 'ABS_y':
        CON = abs(y3_lab) - abs(y4_lab) #abs(y3) - abs(-y3)

    elif FRAME == 'PSE_RAP':
        CON = abs(y3p_lab) - abs(y4p_lab)
        
    elif FRAME == 'OUT_CUT':
        CON1 = abs(y3_lab) - Y_CUT
        CON2 = abs(y4_lab) - Y_CUT

    elif FRAME == 'IN_CUT':
        CON1 = Y_CUT - abs(y4_lab)
        CON2 = Y_CUT - abs(y3_lab)

    #Impose further cuts on y3
    elif FRAME == 'OUT_CUTC':
        CON1 = 0
        CON2 = 0
        if abs(y3_lab) < 2.5:
            CON1 = abs(y3_lab) - Y_CUT
        if abs(y4_lab) < 2.5: 
            CON2 = abs(y4_lab) - Y_CUT
            
    elif FRAME == 'Y_CUT':
        CON1 = 0
        CON2 = 0
        if (y3_lab + y4_lab)/2 > Y_CUT:
            CON1 = y3_lab - y4_lab
            CON2 = y4_lab - y3_lab

    elif FRAME == 'Y_CUTP':
        pT = k3 * sinth * Q * 0.5
        CON1 = 0
        CON2 = 0
        if (y3_lab + y4_lab)/2 > Y_CUT and pT < 20:
            CON1 = y3_lab - y4_lab
            CON2 = y4_lab - y3_lab

    elif FRAME == 'Y_CUTQ':
        CON1 = 0
        CON2 = 0
        if (y3_lab + y4_lab)/2 > Y_CUT and Q > 450:
            CON1 = y3_lab - y4_lab
            CON2 = y4_lab - y3_lab
            
    rf = 0
    rb = 0
    
    if FRAME in FRAME_lB or FRAME in FRAME_lC:
        CON = CON1
    if CON > 0:
        rf += qqV_B(c, hats) - qqV_B(-c, hats) + (1+c**2+4*xm) * \
               (4*np.log((1-c)/(1+c))*np.log(2*w) + qqV_D(c, hats) - qqV_D(-c, hats))
    if FRAME in FRAME_lB or FRAME in FRAME_lC:
        CON = -CON2
    if CON < 0:
        rb += qqV_B(c, hats) - qqV_B(-c, hats) + (1+c**2+4*xm) * \
              (4*np.log((1-c)/(1+c))*np.log(2*w) + qqV_D(c, hats) - qqV_D(-c, hats))

    w_ii = rf - rb
    prefactor = alpha_QCD**3 * 5/108 / hats * pb_convert
    JAC = (1-tau0/x1) * (1 -tau0) * deltath
    w_ii = w_ii * JAC * prefactor * b * fqfqbar_asymmetric(x1, x2) 

    return w_ii

'----=========================================================================----'
'----======================        HARD REGION      ==========================----'
'----=========================================================================----'

def Top_Asym_qqR(x):
    '''Function for the hard real gluon emission'''
    x1 = tau0 + (1 -tau0) * x[0]
    x2 = tau0/x1 +(1-tau0/x1)*x[1]
    hats = x1*x2*s
    Q = math.sqrt(hats)
    xm = M_top**2/hats
    b = math.sqrt(1-4*xm)
    w = E_cut/math.sqrt(hats)

    costh = -1 + deltath*x[2]
    sinth = math.sqrt(1-costh**2)
    phi = 2*math.pi*x[3]
    cosphi = np.cos(phi)
    sinphi = np.sin(phi)

    y12 = 1

    y45_min = w*(1-b)
    y45_max = 1-2*math.sqrt(xm)
    #ADDED CONDITION ELSE LEADS TO NEGATIVE
    if y45_max-y45_min < 0:
        return 0
    y45 = y45_min + (y45_max-y45_min)*x[4]

    E3 = (1-y45)/2
    k3 = math.sqrt(E3**2-xm)
    y13 = E3-k3*costh
    y23 = E3+k3*costh

    y35_min = max((E3-k3)*y45/(1-E3+k3), 2*w-y45)
    y35_max = (E3+k3)*y45/(1-E3-k3)
    y35 = y35_min + (y35_max-y35_min)*x[5]

    E5 = (y35+y45)/2
    if E5< w:
        return 0
    cosa = (y35/(2*E5) -E3) /k3
    sina = math.sqrt(1-cosa**2)
    y34 = 1 - y35 - y45 - 2*xm
    y15 = E5*(1 - sinth*cosphi*sina + costh*cosa)
    y25 = E5*(1 + sinth*cosphi*sina - costh*cosa)
    y14 = 1-y13-y15
    y24 = 1-y23-y25

    E4 = 1-E3-E5
    k4z = -E5*sina*sinth*cosphi + (-k3+E5*cosa)*costh
    y3 = np.log((E3+k3*costh)/(E3-k3*costh)) / 2
    y4 = np.log((E4+k4z)/(E4-k4z)) / 2

    #In lab frame
    y3_lab = y3 + np.log(x1/x2)/2
    y4_lab = y4 + np.log(x1/x2)/2
    #psuedo-rapidity
    k4sq = E4**2 - xm #total momentum k4 squared
    costh4 = k4z/math.sqrt(k4sq)
    y3p = np.log((1 + costh)/(1 - costh))/2
    y4p = np.log((1 + costh4)/(1 - costh4))/2
    y3p_lab = y3p + np.log(x1/x2)/2
    y4p_lab = y4p + np.log(x1/x2)/2
    
    if FRAME == 'ABS_y':
        CON = abs(y3_lab) - abs(y4_lab) # abs(y3) - abs(y4)
        
    elif FRAME == 'PSE_RAP':
        CON = abs(y3p_lab) - abs(y4p_lab)        
        
    elif FRAME == 'OUT_CUT':
        CON1 = abs(y3_lab) - Y_CUT
        CON2 = abs(y4_lab) - Y_CUT
        
    elif FRAME == 'IN_CUT':
        CON1 = Y_CUT - abs(y4_lab)
        CON2 = Y_CUT - abs(y3_lab)
        
    #Impose further cuts on y3
    elif FRAME == 'OUT_CUTC':
        CON1 = 0
        CON2 = 0
        if abs(y3_lab) < 2.5:
            CON1 = abs(y3_lab) - Y_CUT
        if abs(y4_lab) < 2.5: 
            CON2 = abs(y4_lab) - Y_CUT
            
    elif FRAME == 'Y_CUT':
        CON1 = 0
        CON2 = 0
        if (y3_lab + y4_lab)/2 > Y_CUT:
            CON1 = y3_lab - y4_lab
            CON2 = y4_lab - y3_lab

    elif FRAME == 'Y_CUTP':
        pT3 = k3 * sinth * Q * 0.5
        pT4 = math.sqrt(k4sq -k4z**2) * 0.5*Q #transverse momentum of antitop p4
        CON1 = 0
        CON2 = 0
        if (y3_lab + y4_lab)/2 > Y_CUT:
            if pT3 < 20:
                CON1 = y3_lab - y4_lab
            if pT4 < 20:
                CON2 = y4_lab - y3_lab

    elif FRAME == 'Y_CUTQ':
        CON1 = 0
        CON2 = 0
        if (y3_lab + y4_lab)/2 > Y_CUT and Q > 450:
            CON1 = y3_lab - y4_lab
            CON2 = y4_lab - y3_lab
            
    rf = 0
    rb = 0
    
    if FRAME in FRAME_lB or FRAME in FRAME_lC:
        CON = CON1
    if CON > 0:
        rf += ((y13**2 + y14**2 + y23**2 + y24**2 + 2*xm*(y34 + 2*xm + y12)) * y13/y15 + 4*xm*y24) \
               / (y12 * (y34 + 2*xm) * y35) \
               - ((y13**2 + y14**2 + y23**2 + y24**2 + 2*xm*(y34 + 2*xm + y12)) * y23/y25 + 4*xm*y14) \
               / (y12 * (y34 + 2*xm) * y35) \
               - ((y13**2 + y14**2 + y23**2 + y24**2 + 2*xm*(y34 + 2*xm + y12)) * y14/y15 + 4*xm*y23) \
               / (y12 * (y34 + 2*xm) * y45) \
               + ((y13**2 + y14**2 + y23**2 + y24**2 + 2*xm*(y34 + 2*xm + y12)) * y24/y25 + 4*xm*y13) \
               / (y12 * (y34 + 2*xm) * y45)
    if FRAME in FRAME_lB or FRAME in FRAME_lC:
        CON = -CON2
    if CON < 0:
        rb += ((y13**2 + y14**2 + y23**2 + y24**2 + 2*xm*(y34 + 2*xm + y12)) * y13/y15 + 4*xm*y24) \
               / (y12 * (y34 + 2*xm) * y35) \
               - ((y13**2 + y14**2 + y23**2 + y24**2 + 2*xm*(y34 + 2*xm + y12)) * y23/y25 + 4*xm*y14) \
               / (y12 * (y34 + 2*xm) * y35) \
               - ((y13**2 + y14**2 + y23**2 + y24**2 + 2*xm*(y34 + 2*xm + y12)) * y14/y15 + 4*xm*y23) \
               / (y12 * (y34 + 2*xm) * y45) \
               + ((y13**2 + y14**2 + y23**2 + y24**2 + 2*xm*(y34 + 2*xm + y12)) * y24/y25 + 4*xm*y13) \
               / (y12 * (y34 + 2*xm) * y45)
        
    prefactor = alpha_QCD**3 * 5/108 /hats * pb_convert #* 3/8
    JAC = (1-tau0/x1) * (1-tau0) * (y45_max-y45_min) * (y35_max-y35_min) * deltath
    w_ii = (rf - rb) * JAC * prefactor * fqfqbar_asymmetric(x1, x2) 
    return w_ii           

'----=========================================================================----'
'----======================       QUARK-GLUON       ==========================----'
'----=========================================================================----'
def fqfg_asymmetric(x1, x2):
    '''Function for quark-gluon PDFs'''
    mu = M_top
    ff = 0
    quark = [1, 2, 3, 4]
    for ii in range(len(quark)):
        i = quark[ii]
        ff += (( p.xfxQ(i, x1, mu) *  p.xfxQ(21, x2, mu)) - ( p.xfxQ(-i, x1, mu) *  p.xfxQ(21, x2, mu)))
    return ff /(x1 * x2)        
               
def fgfq_asymmetric(x1, x2):
    '''Function for gluon-quark scattering PDFs'''
    mu = M_top
    ff = 0
    quark = [1, 2, 3, 4]
    for ii in range(len(quark)):
        i = quark[ii]
        ff += (( p.xfxQ(21, x1, mu) *  p.xfxQ(i, x2, mu)) - ( p.xfxQ(21, x1, mu) *  p.xfxQ(-i, x2, mu)))
    return ff /(x1 * x2) 

def Top_Asym_qgR(x):
    '''Function foe quark-gluon, gluon-quark scattering'''
    x1 = tau0 + (1 -tau0) * x[0]
    x2 = tau0/x1 +(1-tau0/x1)*x[1]
    hats = x1*x2*s
    Q = math.sqrt(hats)
    xm = M_top**2/hats
    b = math.sqrt(1-4*xm)
    
    costh = -1 + deltath*x[2]
    sinth = math.sqrt(1-costh**2)
    phi = 2*math.pi*x[3]
    cosphi = np.cos(phi)
    sinphi = np.sin(phi)    

    y12 = 1

    y45_min = 0
    y45_max = 1-2*math.sqrt(xm)
    y45 = y45_min + (y45_max-y45_min)*x[4]

    E3 = (1-y45)/2
    k3 = math.sqrt(E3**2-xm)
    y13 = E3-k3*costh
    y23 = E3+k3*costh

    y35_min = (E3-k3)*y45/(1-E3+k3)
    y35_max = (E3+k3)*y45/(1-E3-k3)
    y35 = y35_min + (y35_max-y35_min)*x[5]

    E5 = (y35+y45)/2 
    cosa = (y35/(2*E5) -E3) /k3
    sina = math.sqrt(1-cosa**2)
    y34 = 1 - y35 - y45 - 2*xm
    y15 = E5*(1 - sinth*cosphi*sina + costh*cosa)
    y25 = E5*(1 + sinth*cosphi*sina - costh*cosa)
    y14 = 1-y13-y15
    y24 = 1-y23-y25   
    
    E4 = 1-E3-E5
    k4z = -E5*sina*sinth*cosphi + (-k3+E5*cosa)*costh
    y3qg = np.log((E3+k3*costh)/(E3-k3*costh)) / 2
    y4qg = np.log((E4+k4z)/(E4-k4z)) / 2

    #lab frame
    y3qg_lab = y3qg + np.log(x1/x2)/2
    y4qg_lab = y4qg + np.log(x1/x2)/2
    #psuedo-rapidity
    k4sq = E4**2 - xm #total momentum k4 squared
    costh4 = k4z/math.sqrt(k4sq)
    y3pqg = np.log((1 + costh)/(1 - costh))/2
    y4pqg = np.log((1 + costh4)/(1 - costh4))/2
    y3pqg_lab = y3pqg + np.log(x1/x2)/2
    y4pqg_lab = y4pqg + np.log(x1/x2)/2
    
    if FRAME == 'ABS_y':
        CON = abs(y3qg_lab) - abs(y4qg_lab) #abs(y3qg) - abs(y4qg)

    elif FRAME == 'PSE_RAP':
        CON = abs(y3pqg_lab) - abs(y4pqg_lab)
        
    elif FRAME == 'OUT_CUT':
        CON1 = abs(y3qg_lab) - Y_CUT
        CON2 = abs(y4qg_lab) - Y_CUT
        
    elif FRAME == 'IN_CUT':
        CON1 = Y_CUT - abs(y4qg_lab)
        CON2 = Y_CUT - abs(y3qg_lab)
        
    #Impose further cuts on y3
    elif FRAME == 'OUT_CUTC':
        CON1 = 0
        CON2 = 0
        if abs(y3qg_lab) < 2.5:
            CON1 = abs(y3qg_lab) - Y_CUT
        if abs(y4qg_lab) < 2.5: 
            CON2 = abs(y4qg_lab) - Y_CUT
            
    elif FRAME == 'Y_CUT':
        CON1 = 0
        CON2 = 0
        if (y3qg_lab + y4qg_lab)/2 > Y_CUT:
            CON1 = y3qg_lab - y4qg_lab
            CON2 = y4qg_lab - y3qg_lab

    elif FRAME == 'Y_CUTP':
        pT3 = k3 * sinth * Q * 0.5
        pT4 = math.sqrt(k4sq -k4z**2)* 0.5*Q #transverse momentum of antitop p4
        CON1 = 0
        CON2 = 0
        if (y3qg_lab + y4qg_lab)/2 > Y_CUT:
            if pT3 < 20:
                CON1 = y3qg_lab - y4qg_lab
            if pT4 < 20:
                CON2 = y4qg_lab - y3qg_lab

    elif FRAME == 'Y_CUTQ':
        CON1 = 0
        CON2 = 0
        if (y3qg_lab + y4qg_lab)/2 > Y_CUT and Q > 450:
            CON1 = y3qg_lab - y4qg_lab
            CON2 = y4qg_lab - y3qg_lab
            
    rf1 = 0
    rb1 = 0
    
    if FRAME in FRAME_lB or FRAME in FRAME_lC:
        CON = CON1
    if CON > 0:
        rf1 += ((y13**2 + y14**2 + y35**2 + y45**2 + 2*xm*(y34 + 2*xm -y15)) * (y13/y12 - y35/y25) +
              4*xm*(y45 + y14)) / (y15 * (y34 + 2*xm) * y23) \
              - ((y13**2 + y14**2 + y35**2 + y45**2 + 2*xm*(y34 + 2*xm -y15)) * (y14/y12 - y45/y25) +
                 4*xm*(y35 + y13)) / (y15 *(y34 + 2*xm) * y24)
        
    if FRAME in FRAME_lB or FRAME in FRAME_lC:
        CON = -CON2        
    if CON < 0:
        rb1 += ((y13**2 + y14**2 + y35**2 + y45**2 + 2*xm*(y34 + 2*xm -y15)) * (y13/y12 - y35/y25) +
              4*xm*(y45 + y14)) / (y15 * (y34 + 2*xm) * y23) \
              - ((y13**2 + y14**2 + y35**2 + y45**2 + 2*xm*(y34 + 2*xm -y15)) * (y14/y12 - y45/y25) +
                 4*xm*(y35 + y13)) / (y15 *(y34 + 2*xm) * y24)    

    
    y3gq = np.log((E3 - k3*costh)/(E3 + k3*costh)) / 2
    y4gq = np.log((E4 - k4z)/(E4 + k4z)) / 2
    #lab frame
    y3gq_lab = y3gq + np.log(x1/x2)/2
    y4gq_lab = y4gq + np.log(x1/x2)/2
    #psuedo-rapidity
    y3pgq = np.log((1 - costh)/(1 + costh))/2
    y4pgq = np.log((1 - costh4)/(1 + costh4))/2
    y3pgq_lab = y3pgq + np.log(x1/x2)/2
    y4pgq_lab = y4pgq + np.log(x1/x2)/2
    
    if FRAME == 'ABS_y':
        CON = abs(y3gq_lab) - abs(y4gq_lab) #abs(y3gq) - abs(y4gq)

    elif FRAME == 'PSE_RAP':
        CON = abs(y3pgq_lab) - abs(y4pgq_lab)
        
    elif FRAME == 'OUT_CUT':
        CON1 = abs(y3gq_lab) - Y_CUT
        CON2 = abs(y4gq_lab) - Y_CUT
        
    elif FRAME == 'IN_CUT':
        CON1 = Y_CUT - abs(y4gq_lab)
        CON2 = Y_CUT - abs(y3gq_lab)
        
    #Impose further cuts on y3
    elif FRAME == 'OUT_CUTC':
        CON1 = 0
        CON2 = 0
        if abs(y3gq_lab) < 2.5:
            CON1 = abs(y3gq_lab) - Y_CUT
        if abs(y4gq_lab) < 2.5: 
            CON2 = abs(y4gq_lab) - Y_CUT
            
    elif FRAME == 'Y_CUT':
        CON1 = 0
        CON2 = 0
        if (y3gq_lab + y4gq_lab)/2 > Y_CUT:
            CON1 = y3gq_lab - y4gq_lab
            CON2 = y4gq_lab - y3gq_lab

    elif FRAME == 'Y_CUTP':
        pT3 = k3 * sinth * Q * 0.5
        pT4 = math.sqrt(k4sq -k4z**2) * 0.5 *Q #transverse momentum of antitop p4
        CON1 = 0
        CON2 = 0
        if (y3gq_lab + y4gq_lab)/2 > Y_CUT:
            if pT3 < 20:
                CON1 = y3gq_lab - y4gq_lab
            if pT4 < 20:
                CON2 = y4gq_lab - y3gq_lab

    elif FRAME == 'Y_CUTQ':
        CON1 = 0
        CON2 = 0
        if (y3gq_lab + y4gq_lab)/2 > Y_CUT and Q > 450:
            CON1 = y3gq_lab - y4gq_lab
            CON2 = y4gq_lab - y3gq_lab
            
    rf2 = 0
    rb2 = 0
    
    if FRAME in FRAME_lB or FRAME in FRAME_lC:
        CON = CON1 
    if CON > 0:
        rf2 += ((y13**2 + y14**2 + y35**2 + y45**2 + 2*xm*(y34 + 2*xm -y15)) * (y13/y12 - y35/y25) +
              4*xm*(y45 + y14)) / (y15 * (y34 + 2*xm) * y23) \
              - ((y13**2 + y14**2 + y35**2 + y45**2 + 2*xm*(y34 + 2*xm -y15)) * (y14/y12 - y45/y25) +
                 4*xm*(y35 + y13)) / (y15 *(y34 + 2*xm) * y24)
        
    if FRAME in FRAME_lB or FRAME in FRAME_lC:
        CON = -CON2
    if CON < 0:
        rb2 += ((y13**2 + y14**2 + y35**2 + y45**2 + 2*xm*(y34 + 2*xm -y15)) * (y13/y12 - y35/y25) +
              4*xm*(y45 + y14)) / (y15 * (y34 + 2*xm) * y23) \
              - ((y13**2 + y14**2 + y35**2 + y45**2 + 2*xm*(y34 + 2*xm -y15)) * (y14/y12 - y45/y25) +
                 4*xm*(y35 + y13)) / (y15 *(y34 + 2*xm) * y24)
        
    prefactor = alpha_QCD**3 * 5/108 / hats * pb_convert * 3/8
    JAC = (1-tau0/x1) * (1-tau0) * (y45_max-y45_min) * (y35_max-y35_min) * deltath
    w_ii = (rf1-rb1) * JAC * prefactor * fqfg_asymmetric(x1, x2)
    w_ii +=  (rf2-rb2) * JAC * prefactor * fgfq_asymmetric(x1, x2) 

    return w_ii
'----=========================================================================----'
'----=======================     ALPHA QCD COUPLING   ========================----'
'----=========================================================================----'


'----=========================================================================----'
'----=========================================================================----'
'----=======================     INTEGRATE ASYMM     =========================----'
'----=========================================================================----'
'----=========================================================================----'

 #make directory path
date = time.strftime("%m-%d-%H:%M")
cwd = os.getcwd()
base = os.path.splitext(os.path.basename(__file__))[0]
print date


 #make generic dict
dictt = {} #generic dictionary with mu = M_top

 #PDF sets
##PDF0 = 'cteq6l1'
##CTEQ = 'cteq66'
##MSTW_LO = 'MSTW2008lo90cl'
MSTW_NLO = 'MSTW2008nlo90cl'
##NNPDF = 'NNPDF21_lo_as_0119_100'
##NNPDF_STAR = 'NNPDF21_lostar_as_0119_100'
PDF = MSTW_NLO
 # initializes PDF member object (for protons)
p = lhapdf.mkPDF(PDF, 0)

 #Frame = ['ABS_y', 'OUT_CUT', 'IN_CUT', 'OUT_CUTC', 'Y_CUT',  'Y_CUTP', 'Y_CUTQ','PSE_RAP']
FRAME_lA = ['ABS_y','PSE_RAP']
FRAME_lB = ['OUT_CUT', 'IN_CUT', 'OUT_CUTC']
FRAME_lC = ['Y_CUT',  'Y_CUTP', 'Y_CUTQ']

NEVAL = 1000000
NEV_CON = 4
print ('NEVAL = %s' % NEVAL)

FRAME = 'OUT_CUTC'
print FRAME
mu_list = [M_top/2] ##M_top, M_top/2, M_top *2
name_list = ['M_half'] #'M_topp', 'M_half', 'M_doub'


if FRAME in FRAME_lA:
    # PP COM energy in GeV
    E_list = [7000, 8000, 10000, 12000, 14000]
    name_ABS = ['7TeV', '8TeV', '10TeV', '12TeV', '14TeV']
    
if FRAME in FRAME_lB:
##    CUT_list = [0.25, 0.35, 0.5, 0.65, 0.75, 0.85, 1.0, 1.15, 1.25, 1.40, 1.5, 1.60, 1.75, 1.80, 2.0, 2.2, 2.4, 2.5]
    CUT_list = [0.05, 0.25, 0.50 ,0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50] #  
    # PP COM energy in GeV
    ECM = 7000 #set energy
    s = ECM**2
    
if FRAME in FRAME_lC:
##    CUT_list = [0.25, 0.5, 0.70, 1.0, 1.40, 1.70, 2.0, 2.2, 2.5, 2.55, 2.6, 2.65, 2.7, 2.9]
##    name_CUT = ['0.25', '0.5', '0.70', '1.0', '1.40', '1.70', '2.0', '2.2', '2.5', '2.55', '2.6', '2.65', '2.7', '2.9']
    CUT_list = [1.50, 1.55, 1.60, 1.65, 1.70, 1.80, 1.95, 2.0, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.40, 2.45, 2.50]
    # PP COM energy in GeV
    ECM = 7000 #set energy
    s = ECM**2

'----=========================================================================----'
'----=======================     SAVE FILE SETUP     =========================----'
'----=========================================================================----'
if NEVAL > NEV_CON:
    Save_path1 = cwd+ '/' + base[0:3] +'/SAV_'+base+'_'+FRAME+'/%s/' % (date)
    Save_path = Save_path1
    if not os.path.exists(os.path.dirname(Save_path)):
        try:
            os.makedirs(os.path.dirname(Save_path))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

     #save text
    with open(os.path.join(Save_path1, 'details.txt'), 'w') as myfile:
        myfile.write('\n Array variables \n')


for kk in range(len(mu_list)):
    mu = mu_list[kk]
    alpha_QCD = p.alphasQ(mu)
    print ('mu (GeV) = %s alpha_QCD = %s' % (mu, alpha_QCD))
    print date
    '----=========================================================================----'
    '----=======================     SAVE FILE SETUP     =========================----'
    '----=========================================================================----'
    if NEVAL > NEV_CON:
        Save_path1 = cwd+ '/' + base[0:3] +'/SAV_'+base+'_'+FRAME+'/%s/' % (date)
        Save_path = Save_path1
        if not os.path.exists(os.path.dirname(Save_path)):
            try:
                os.makedirs(os.path.dirname(Save_path))
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise

         #save text
        with open(os.path.join(Save_path1, 'details.txt'), 'a') as myfile:
            myfile.write('\n====================================================================================\n')
            myfile.write('====================================================================================\n')
            myfile.write('\n|| Frame = %s || mu_f = %s || alphaQCD = %s || NEVAL = %s ||\n' % (FRAME, mu, alpha_QCD, NEVAL))
                    
 #Choose number of iterations
    if FRAME in FRAME_lA:
        iterations = len(E_list)
    elif FRAME in FRAME_lB or FRAME in FRAME_lC:
        iterations = len(CUT_list)

    for i in range(iterations):
        print '----=========================================================================----'
        print '----=========================================================================----'
        print time.strftime("%H:%M")
        t0 = time.strftime("%H:%M")

        if FRAME in FRAME_lB or FRAME in FRAME_lC:
            Y_CUT = CUT_list[i]
            name = '%.2f' % Y_CUT
            
        if FRAME in FRAME_lA:
            ECM = E_list[i]
            s = ECM**2
            name = name_ABS[i]
            
        print "hadron com energy:", ECM, "GeV"
        tau0 = 4*M_top**2/s


        print ('mu (GeV) = %s alpha_QCD = %s and Frame: %s' % (mu, alpha_QCD, FRAME))
        with open(os.path.join(Save_path1, 'details.txt'), 'a') as myfile:
            myfile.write('\n =================================================')            
            if FRAME in FRAME_lB or FRAME in FRAME_lC:                
                myfile.write('\n || Y cut = %s ||\n' % Y_CUT)
            elif FRAME in FRAME_lA:
                myfile.write('\n || ECM = %s ||\n' % ECM)
        mem()
        '----=========================================================================----'
        '----========================     qq CALCULATION    ==========================----'
        '----=========================================================================----'
        
        print time.strftime("%H:%M")
        print 'Integrating qq no. %s for mu = %s' % (i+1, mu)
        integ_LO_qq = vegas.Integrator(3 * [[0,1]])
        result_LO_qq = integ_LO_qq(Top_LO_qq, nitn=10, neval= NEVAL)
        print(result_LO_qq.summary())
        print('sigma CS(pb) = %s    Q = %.2f' % (result_LO_qq, result_LO_qq.Q))

        reg = 'qqbar'
        #create lists for x, wgt in integ
        x0 = []
        x1 = []
        x2 = []
        wghts = []

        summed_res = 0
        jj = 0
        Neve = sum(1 for i in integ_LO_qq.random())
        for x, wgt in integ_LO_qq.random():
            summed_res += wgt * Top_LO_qq(x)
            sys.stdout.write("progress: %d%%   \r" % (float(jj)*100./(Neve)) )
            sys.stdout.flush()
##            x0.append(x[0])
##            x1.append(x[1])
##            x2.append(x[2])
            wghts.append(Top_LO_qq(x) * wgt)
            jj+=1

        print ('LO qq CS(pb) = %s' % summed_res)
        #create dictionary
        qq_bar = {'x0': x0, 'x1': x1, 'x2': x2, 'wghts': wghts}

        #append into appropriate dictionary for mu_f

        if NEVAL > NEV_CON:
            pickle_out = open(Save_path + name_list[kk]+'_' +FRAME +'_' + name +'_'+ reg+'.pickle', 'wb')
            pickle.dump(qq_bar, pickle_out)
            pickle_out.close()

        if NEVAL > NEV_CON:
            with open(os.path.join(Save_path1, 'details.txt'), 'a') as myfile:
                myfile.write(' || LO_qq CS(pb) = %s \n' % (summed_res))

        mem()
        print 'clearing array ...'
        qq_bar = {}
        gc.collect()
        x0 = []
        x1 = []
        x2 = []
        wghts = []
        mem()


        '----=========================================================================----'
        '----========================     gg CALCULATION    ==========================----'
        '----=========================================================================----'
        
        print time.strftime("%H:%M")
        print ('Integrating gg no. %s with mu = %s' % (i+1, mu))
        integ_LO_gg = vegas.Integrator(3 * [[0,1]])
        result_LO_gg = integ_LO_gg(Top_LO_gg, nitn=10, neval= NEVAL)
        print(result_LO_gg.summary())
        print('sigma CS(pb) = %s    Q = %.2f' % (result_LO_gg, result_LO_gg.Q))

        reg = 'gg'
        #create lists for x, wgt in integ
        x0 = []
        x1 = []
        x2 = []
        wghts = []

        summed_res = 0
        jj = 0
        Neve = sum(1 for i in integ_LO_gg.random())
        for x, wgt in integ_LO_gg.random():
            summed_res += wgt * Top_LO_gg(x)
            sys.stdout.write("progress: %d%%   \r" % (float(jj)*100./(Neve)) )
            sys.stdout.flush()
##            x0.append(x[0])
##            x1.append(x[1])
##            x2.append(x[2])
            wghts.append(Top_LO_gg(x) * wgt)
            jj+=1
            
        print ('LO gg CS(pb) = %s' % summed_res)
        #create dictionary    
        gg = {'x0': x0, 'x1': x1, 'x2': x2, 'wghts': wghts}
        #append into appropriate dictionary for mu_f
        if NEVAL > NEV_CON:
            pickle_out = open(Save_path + name_list[kk]+'_' +FRAME +'_' + name +'_'+ reg+'.pickle', 'wb')
            pickle.dump(gg, pickle_out)
            pickle_out.close()

        if NEVAL > NEV_CON:
            with open(os.path.join(Save_path1, 'details.txt'), 'a') as myfile:
                myfile.write(' || LO_gg CS(pb) = %s \n' % (summed_res))

        mem()
        print 'clearing array ...'
        gg = {}
##        x0 = []
##        x1 = []
##        x2 = []
        wghts = []
        gc.collect()
        mem()

        '----=========================================================================----'
        '----========================       SOFT REGION     ==========================----'
        '----=========================================================================----'
        
        print time.strftime("%H:%M")
        print 'Integrating soft no. %s for mu = %s' % (i+1, mu)
        integ_soft = vegas.Integrator(3 * [[0,1]])
        result_soft = integ_soft(Top_Asym_qqV, nitn=10, neval= NEVAL) ##
        print(result_soft.summary())
        print('sigma CS(pb) = %s    Q = %.2f' % (result_soft, result_soft.Q))

        reg = 'soft'
        #create lists for x, wgt in integ
        x0 = []
        x1 = []
        x2 = []
        wghts = []

        summed_res = 0
        jj = 0
        Neve = sum(1 for i in integ_soft.random())
        for x, wgt in integ_soft.random():
            summed_res += Top_Asym_qqV(x) * wgt
            sys.stdout.write("progress: %d%%   \r" % (float(jj)*100./(Neve)) )
            sys.stdout.flush()
##            x0.append(x[0])
##            x1.append(x[1])
##            x2.append(x[2])
            wghts.append(Top_Asym_qqV(x) * wgt)
            jj+=1
            
        print ('soft CS(pb) = %s' % summed_res)
        #create dictionary
        soft = {'x0': x0, 'x1': x1, 'x2': x2, 'wghts': wghts}
        #append into appropriate dictionary for mu_f
        if NEVAL > NEV_CON:
            pickle_out = open(Save_path + name_list[kk]+'_' +FRAME +'_' + name +'_'+ reg+'.pickle', 'wb')
            pickle.dump(soft, pickle_out)
            pickle_out.close()

        if NEVAL > NEV_CON:
            with open(os.path.join(Save_path1, 'details.txt'), 'a') as myfile:
                myfile.write(' || soft CS(pb) = %s \n' % (summed_res))

        mem()
        print 'clearing array ...'
        soft = {}
        x0 = []
        x1 = []
        x2 = []
        wghts = []    
        gc.collect()
        mem()

        '----=========================================================================----'
        '----========================       HARD REGION     ==========================----'
        '----=========================================================================----'
        
        print time.strftime("%H:%M")
        print 'Integrating hard no. %s for mu = %s' % (i+1, mu)
        integ_hard = vegas.Integrator(6 * [[0,1]])
        result_hard = integ_hard(Top_Asym_qqR, nitn=10, neval= NEVAL) ##
        print(result_hard.summary())
        print('sigma CS(pb) = %s    Q = %.2f' % (result_hard, result_hard.Q))

        reg = 'hard'
        #create lists for x, wgt in integ
        x0 = []
        x1 = []
        x2 = []
        x3 = []
        x4 = []
        x5 = []
        wghts = []

        summed_res = 0
        jj = 0
        Neve = sum(1 for i in integ_hard.random())
        for x, wgt in integ_hard.random():
            summed_res += Top_Asym_qqR(x) * wgt
            sys.stdout.write("progress: %d%%   \r" % (float(jj)*100./(Neve)) )
            sys.stdout.flush()
##            x0.append(x[0])
##            x1.append(x[1])
##            x2.append(x[2])
##            x3.append(x[3])
##            x4.append(x[4])
##            x5.append(x[5])
            wghts.append(Top_Asym_qqR(x) * wgt)
            jj+=1

        print ('hard CS(pb) = %s' % summed_res)
        #create dictionary
        hard = {'x0': x0, 'x1': x1, 'x2': x2, 'x3': x3, 'x4': x4, 'x5': x5, 'wghts': wghts}

        #append into appropriate dictionary for mu_f
        if NEVAL > NEV_CON:
            pickle_out = open(Save_path + name_list[kk]+'_' +FRAME +'_' + name +'_'+ reg+'.pickle', 'wb')
            pickle.dump(hard, pickle_out)
            pickle_out.close()

        if NEVAL > NEV_CON:
            with open(os.path.join(Save_path1, 'details.txt'), 'a') as myfile:
                myfile.write(' || hard CS(pb) = %s \n' % (summed_res))

        mem()
        print 'clearing array ...'
        hard = {}
        x0 = []
        x1 = []
        x2 = []
        x3 = []
        x4 = []
        x5 = []
        wghts = []    
        gc.collect()
        mem()

        '----=========================================================================----'
        '----========================        QG REGION      ==========================----'
        '----=========================================================================----'

        print time.strftime("%H:%M")
        print 'Integrating quark-gluon no. %s for mu = %s' % (i+1, mu)
        integ_qg = vegas.Integrator(6 * [[0,1]])
        result_qg = integ_qg(Top_Asym_qgR, nitn=10, neval= NEVAL) ##
        print(result_qg.summary())
        print('sigma CS(pb) = %s    Q = %.2f' % (result_qg, result_qg.Q))

        reg = 'qg'
        #create lists for x, wgt in integ
        x0 = []
        x1 = []
        x2 = []
        x3 = []
        x4 = []
        x5 = []
        wghts = []

        summed_res = 0
        jj = 0
        Neve = sum(1 for i in integ_qg.random())
        for x, wgt in integ_qg.random():
            summed_res += Top_Asym_qgR(x) * wgt
            sys.stdout.write("progress: %d%%   \r" % (float(jj)*100./(Neve)) )
            sys.stdout.flush()
##            x0.append(x[0])
##            x1.append(x[1])
##            x2.append(x[2])
##            x3.append(x[3])
##            x4.append(x[4])
##            x5.append(x[5])
            wghts.append(Top_Asym_qgR(x) * wgt)
            jj+=1

        print ('qg CS(pb) = %s' % summed_res)
        #create dictionary
        qg = {'x0': x0, 'x1': x1, 'x2': x2, 'x3': x3, 'x4': x4, 'x5': x5, 'wghts': wghts}

            #append into appropriate dictionary for mu_f
        if NEVAL > NEV_CON:
            pickle_out = open(Save_path + name_list[kk]+'_' +FRAME +'_' + name +'_'+ reg+'.pickle', 'wb')
            pickle.dump(qg, pickle_out)
            pickle_out.close()

        t1 = time.strftime("%H:%M")

        if NEVAL > NEV_CON:
            with open(os.path.join(Save_path1, 'details.txt'), 'a') as myfile:
                myfile.write(' || quark-gluon CS(pb) = %s \n' % (summed_res))
                myfile.write(' || time elapsed: %s - %s\n' % (t0, t1))

        mem()
        print 'clearing array ...'
        qg = {}
        x0 = []
        x1 = []
        x2 = []
        x3 = []
        x4 = []
        x5 = []
        wghts = []
        gc.collect()
        mem()

        
Binning_file = '/home/chiho/Documents/Chiho Documents/Physics University Durham/4th year/Final Project/Computing/Top_Quark/Binning_tools/Bin_Obs_LHC.py'
Bin_dir = Save_path = cwd+ '/' +base[0:3]+'/SAV_'+base+'_'+FRAME+'/%s/Bin_Obs_LHC.py' % date
if NEVAL > NEV_CON:
    shutil.copy(Binning_file, Bin_dir)

print date[-5:] + ' - ' + time.strftime("%H:%M")  


