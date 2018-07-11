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
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import sys
from scipy.special import spence as Li
import time
import errno

sys.path.append('/home/chiho/lhapdf/lib/python2.7/site-packages')

import lhapdf
# initializes PDF member object (for protons)
p = lhapdf.mkPDF("cteq6l1", 0)

## some parameters
pb_convert = 3.894E8 # conversion factor GeV^-2 -> pb
MZ = 91.188    # Z boson mass
GAMMAZ = 2.4414 # Z boson width
alpha = 1/132.507 # alpha QED
Gfermi = 1.16639E-5 # fermi constant
sw2 = 0.222246 # sin^2(weinberg angle)

#Define variables in Kuhn paper
M_top = 173.1 # Mass of top quark`
alpha_QCD = 0.109 # alpha QCD
##E_cut = 0.001

# PP COM energy in GeV
ECM = 1960
s = ECM**2
print "hadron com energy:", ECM, "GeV"

# we also need the "ranges" (i.e. x2 - x1)
deltath = 2
tau0 = 4*M_top**2/s

print '\n'
print '----====================================================----'
print 'Exercise 2 (Advanced Scientific Computing Workshop, July 2014)'
print 'pp_bar --> g --> t t_bar'
print '----====================================================----'
print '\n'

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
        ff += (p.xfxQ(i, x1, mu) *  p.xfxQ(i, x2, mu)) + (p.xfxQ(-i, x1, mu) *  p.xfxQ(-i, x2, mu))
    return ff /(x1 * x2)

def fgfg_symmetric(x1, x2):
    '''Function for the gg symmetric PDFs'''
    mu = M_top
    ff = 0
    ff += (p.xfxQ(21, x1, mu) *  p.xfxQ(21, x2, mu))  # + (p.xfxQ(-i, x1, mu) *  p.xfxQ(-i, x2, mu))
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
    
    #random costh value
    costh = -1 + deltath*x[2]
    beta = math.sqrt(1 - (4 * xm))
    c = beta * costh

    JAC = (1 -tau0) * (1-tau0/x1) * deltath
    PREFAC = alpha_QCD**2 * math.pi * beta * pb_convert /(9 * hats) 
    return  PREFAC * JAC *(1 + c**2 + (4 * xm)) * fqfqbar_symmetric(x1, x2)

def Top_LO_gg(x):
    x1 = tau0 + (1 -tau0) * x[0]
    x2 = tau0/x1 +(1-tau0/x1)*x[1]
    hats = x1*x2*s
    xm = M_top **2 / hats
    if hats < 4*M_top**2:
        print 'hats < 4M_top'
    #random costh value
    costh = -1 + deltath*x[2]
    beta = math.sqrt(1 - (4 * xm))
    c = beta * costh

    JAC = (1 -tau0) * (1-tau0/x1) * deltath
    PREFAC = alpha_QCD**2 * math.pi * beta * pb_convert /(2 * hats)
    w_ii = (1/(3 * (1-c**2)) - 3/16) * (1 + c**2 + 8*xm - (32*xm**2)/(1-c**2))
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
        ff += (p.xfxQ(i, x1, mu) *  p.xfxQ(i, x2, mu)) - (p.xfxQ(-i, x1, mu) *  p.xfxQ(-i, x2, mu))
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
    b = math.sqrt(1-4*xm)
    w = E_cut/math.sqrt(hats)

    costh = -1 + deltath*x[2]
    c = b*costh

    E3 = 0.5
    k3 = math.sqrt(E3**2 -xm)
    y3 = np.log((E3 + k3*costh)/(E3 - k3*costh))/2
    y4 = np.log((E3 - k3*costh)/(E3 + k3*costh))/2
    y3_lab = y3 + np.log(x1/x2)/2
    
    if y3_lab < 0:
        return 0

    w_ii = qqV_B(c, hats) - qqV_B(-c, hats) + (1+c**2+4*xm) * \
           (4*np.log((1-c)/(1+c))*np.log(2*w) + qqV_D(c, hats) - qqV_D(-c, hats))
    prefactor = alpha_QCD**3 * 5/108 / hats * pb_convert
    JAC = (1-tau0/x1) * (1 -tau0) * deltath
    w_ii = w_ii * JAC * prefactor * b * fqfqbar_asymmetric(x1, x2) *2 
    return w_ii

'----=========================================================================----'
'----======================        HARD REGION      ==========================----'
'----=========================================================================----'

def Top_Asym_qqR(x):
    '''Function for the hard real gluon emission'''
    x1 = tau0 + (1 -tau0) * x[0]
    x2 = tau0/x1 +(1-tau0/x1)*x[1]
    hats = x1*x2*s
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
    y3_lab = y3 + np.log(x1/x2)/2
    
    y_cond =  y3-y4
    if y3_lab < 0:
        return 0
    
    w_ii = ((y13**2 + y14**2 + y23**2 + y24**2 + 2*xm*(y34 + 2*xm + y12)) * y13/y15 + 4*xm*y24) \
           / (y12 * (y34 + 2*xm) * y35) \
           - ((y13**2 + y14**2 + y23**2 + y24**2 + 2*xm*(y34 + 2*xm + y12)) * y23/y25 + 4*xm*y14) \
           / (y12 * (y34 + 2*xm) * y35) \
           - ((y13**2 + y14**2 + y23**2 + y24**2 + 2*xm*(y34 + 2*xm + y12)) * y14/y15 + 4*xm*y23) \
           / (y12 * (y34 + 2*xm) * y45) \
           + ((y13**2 + y14**2 + y23**2 + y24**2 + 2*xm*(y34 + 2*xm + y12)) * y24/y25 + 4*xm*y13) \
           / (y12 * (y34 + 2*xm) * y45)
    prefactor = alpha_QCD**3 * 5/108 /hats * pb_convert
    JAC = (1-tau0/x1) * (1-tau0) * (y45_max-y45_min) * (y35_max-y35_min) * deltath
    w_ii = w_ii * JAC * prefactor * fqfqbar_asymmetric(x1, x2) *2
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
        ff += (( p.xfxQ(21, x1, mu) *  p.xfxQ(-i, x2, mu)) - ( p.xfxQ(21, x1, mu) *  p.xfxQ(i, x2, mu)))
    return ff /(x1 * x2) 

def Top_Asym_qgR(x):
    '''Function foe quark-gluon, gluon-quark scattering'''
    x1 = tau0 + (1 -tau0) * x[0]
    x2 = tau0/x1 +(1-tau0/x1)*x[1]
    hats = x1*x2*s
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
    y3qg_lab = y3qg + np.log(x1/x2)/2

    t1 = 0
    t2 = 0

    y3qg_delta = y3qg-y4qg
    if y3qg_lab > 0:
        t1 += ((y13**2 + y14**2 + y35**2 + y45**2 + 2*xm*(y34 + 2*xm -y15)) * (y13/y12 - y35/y25) +
              4*xm*(y45 + y14)) / (y15 * (y34 + 2*xm) * y23) \
              - ((y13**2 + y14**2 + y35**2 + y45**2 + 2*xm*(y34 + 2*xm -y15)) * (y14/y12 - y45/y25) +
                 4*xm*(y35 + y13)) / (y15 *(y34 + 2*xm) * y24)

    y3gq = np.log((E3 - k3*costh)/(E3 + k3*costh)) / 2
    y4gq = np.log((E4 - k4z)/(E4 + k4z)) / 2
    y3gq_lab = y3gq + np.log(x1/x2)/2
    
    y3gq_delta = y3gq-y4gq
    if y3gq_lab > 0:
        t2 += ((y13**2 + y14**2 + y35**2 + y45**2 + 2*xm*(y34 + 2*xm -y15)) * (y13/y12 - y35/y25) +
              4*xm*(y45 + y14)) / (y15 * (y34 + 2*xm) * y23) \
              - ((y13**2 + y14**2 + y35**2 + y45**2 + 2*xm*(y34 + 2*xm -y15)) * (y14/y12 - y45/y25) +
                 4*xm*(y35 + y13)) / (y15 *(y34 + 2*xm) * y24)

    prefactor = alpha_QCD**3 * 5/108 / hats * pb_convert
    JAC = (1-tau0/x1) * (1-tau0) * (y45_max-y45_min) * (y35_max-y35_min) * deltath
    w_ii = t1 * JAC * prefactor * fqfg_asymmetric(x1, x2) *2
    w_ii = w_ii + JAC * t2 * prefactor * fgfq_asymmetric(x1, x2) *2 

    return w_ii

'----=========================================================================----'
'----=========================================================================----'
'----================      Virtual + SOFT & REAL SUM      ====================----'
'----=========================================================================----'
'----=========================================================================----'

#make directory path
date = time.strftime("%m-%d-%H:%M")
cwd = os.getcwd()
base = os.path.splitext(os.path.basename(__file__))[0]
print date
NEVAL = 500000
NEV_CON = 9000
print ('NEVAL = %s' % NEVAL)
if NEVAL > NEV_CON:
    Save_path = cwd+'/SAV_DAT/'
    if not os.path.exists(os.path.dirname(Save_path)):
        try:
            os.makedirs(os.path.dirname(Save_path))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

print time.strftime("%m-%d-%H:%M")
E_cut_list =[0.00001, 0.0001, 0.001, 0.01, 0.1]
E_cut_list = [i * M_top for i in E_cut_list]
if NEVAL > NEV_CON:
    SaveStr = Save_path +'/E_cut'
    np.save(SaveStr, E_cut_list)

'----=========================================================================----'
'----========================       SOFT REGION     ==========================----'
'----=========================================================================----'

print E_cut_list, 'E_cut list'
soft_CS = []
for i in range(len(E_cut_list)):
    ''' Integrating Soft region for a range of E_cuts'''
    print time.strftime("%H:%M")
    E_cut = E_cut_list[i]
    print ('Integrating Soft region no. %s with E_cut = %s' % (i+1, E_cut))
    integ_soft = vegas.Integrator(3 * [[0,1]])
    result_soft = integ_soft(Top_Asym_qqV, nitn=10, neval= NEVAL)
    print(result_soft.summary())
    print('sigma CS(pb) = %s    Q = %.2f' % (result_soft, result_soft.Q))    
    summed_res = 0
    Neve = sum(1 for i in integ_soft.random())
    jj = 0
    for x, wgt in integ_soft.random():
        sys.stdout.write("progress: %d%%   \r" % (float(jj)*100./(Neve)) )
        sys.stdout.flush()
        summed_res += wgt * Top_Asym_qqV(x)
        jj += 1
    print summed_res, 'summed res'    
    soft_CS.append(summed_res)
    if NEVAL > NEV_CON:
        with open(os.path.join(Save_path, 'soft_val.txt'), 'a') as myfile:
            myfile.write('\n || Soft CS(pb) = %s || E cut = %s || NEVAL = %s ||' % (summed_res, E_cut, NEVAL))

if NEVAL > NEV_CON:
    SaveStr = Save_path + '/soft_CS'        
    np.save(SaveStr, soft_CS)

'----=========================================================================----'
'----========================       HARD REGION     ==========================----'
'----=========================================================================----'

hard_CS = []
for i in range(len(E_cut_list)):
    ''' Integrating Hard region for a range of E_cuts'''
    print time.strftime("%H:%M")
    E_cut = E_cut_list[i]
    print ('Integrating Hard region no. %s with E_cut = %s' % (i+1, E_cut))
    integ_hard = vegas.Integrator(6 * [[0,1]])
    result_hard = integ_hard(Top_Asym_qqR, nitn=10, neval= NEVAL)
    print(result_hard.summary())
    print('sigma CS(pb) = %s    Q = %.2f' % (result_hard, result_hard.Q))    
    summed_res = 0
    Neve = sum(1 for i in integ_hard.random())
    jj = 0
    for x, wgt in integ_hard.random():
        sys.stdout.write("progress: %d%%   \r" % (float(jj)*100./(Neve)) )
        sys.stdout.flush()
        summed_res += wgt * Top_Asym_qqR(x)
        jj += 1
    print summed_res, 'summed res'    
    hard_CS.append(summed_res)
    if NEVAL > NEV_CON:
        with open(os.path.join(Save_path, 'hard_val.txt'), 'a') as myfile:
            myfile.write('\n || Hard CS(pb) = %s || E cut = %s || NEVAL = %s ||' % (summed_res, E_cut, NEVAL))

if NEVAL > NEV_CON:
    SaveStr = Save_path + '/hard_CS'
    np.save(SaveStr, hard_CS)

tot_hs_CS = [x + y for x, y in zip(soft_CS, hard_CS)]
if NEVAL > NEV_CON:
    SaveStr = Save_path + '/tot_hs_CS'
    np.save(SaveStr, tot_hs_CS)
print time.strftime("%H:%M")

'----=========================================================================----'
'----========================        QG REGION      ==========================----'
'----=========================================================================----'
    
''' Integrating quark-gluon channel'''
print 'Integrating quark-gluon scattering'
integ_qg = vegas.Integrator(6 * [[0,1]])
result_qg = integ_qg(Top_Asym_qgR, nitn=10, neval= NEVAL)
print(result_qg.summary())
print('sigma CS(pb) = %s    Q = %.2f' % (result_qg, result_qg.Q))    
qg_CS = 0
Neve_qg = sum(1 for i in integ_qg.random())
jj = 0
for x, wgt in integ_qg.random():
    sys.stdout.write("progress: %d%%   \r" % (float(jj)*100./(Neve_qg)) )
    sys.stdout.flush()
    qg_CS += wgt * Top_Asym_qgR(x)
    jj += 1
print qg_CS, 'qg summed res'

if NEVAL > NEV_CON:
    with open(os.path.join(Save_path, 'qg.txt'), 'a') as myfile:
        myfile.write('\n || Hard CS(pb) = %s || E cut = %s || NEVAL = %s ||' % (summed_res, E_cut, NEVAL))

tot_hsq_CS = [x + y + qg_CS for x, y in zip(soft_CS, hard_CS)]

if NEVAL > NEV_CON:
    SaveStr = Save_path + '/tot_hsq_CS'
    np.save(SaveStr, tot_hsq_CS)
print time.strftime("%H:%M")
    
##plt.figure()
##plt.semilogx(E_cut_list, soft_CS, label = 'soft', color='blue')
##plt.semilogx(E_cut_list, hard_CS, label = 'hard', color='green')
##plt.semilogx(E_cut_list, tot_hsq_CS, linestyle='--', label = 'hard', color='red')
##plt.xlabel(r'$E_{cut}$')
##plt.ylabel(r'$\sigma^{Tot}_{A}$')
##
##plt.tight_layout()
##plt.show()

