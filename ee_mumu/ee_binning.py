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
import glob
from matplotlib import gridspec
from matplotlib.ticker import AutoMinorLocator

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

'----=========================================================================----'
'----=======================     HISTOGRAM BINNING    ========================----'
'----=========================================================================----'

## #what follows depends on matplotlib/numpy
## #(taken from http://matplotlib.org/examples/api/histogram_path_demo.html)
def generate_histo(array, name, xlabel, array_range = 0):
  
    fig, ax = plt.subplots()
    if array_range is 0:
        array_range = (np.min(array), np.max(array))

    n, bins = np.histogram(array, 50, array_range)
    
    # get the corners of the rectangles for the histogram
    left = np.array(bins[:-1])
    right = np.array(bins[1:])
    bottom = np.zeros(len(left))
    top = (bottom + n)

    # function to build a compound path
    XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T
    
    # get the Path object
    barpath = path.Path.make_compound_path_from_polys(XY)
    
    # make a patch out of it
    patch = patches.PathPatch(barpath, facecolor=colour, edgecolor='black', alpha=1)
    ax.add_patch(patch)
    
    #draw legend
    legend_patch = [patches.Patch(color=colour, label = "{:s}".format(r'$e^{-}e^{+}$'r'$\rightarrow$ Z'r'/$\gamma$'r'$\rightarrow$'r'$\mu^{-}\mu^{-}$'))]
    plt.legend(handles=legend_patch, bbox_to_anchor=(0.97, 0.99))
    
    # update the view limits
    ax.set_xlim(left[0], right[-1])
    ax.set_ylim(bottom.min(), top.max() + 0.08*top.max())
    plt.xlabel(xlabel)
    plt.ylabel('Number of events')
    plt.savefig(Save_path + name + '.pdf')
##    plt.show()

def generate_normalized_histo(array, wghts, name, xlabel, ylabel, array_range = 0):
    print ylabel
    fig, ax = plt.subplots()
    if array_range is 0:
        array_range = (np.min(array), np.max(array))
    n, bins = np.histogram(array, 50, weights=wghts)
    
    # get the corners of the rectangles for the histogram                                                                
    left = np.array(bins[:-1])
    right = np.array(bins[1:])
    bottom = np.zeros(len(left))
    top = (bottom + n) / ((right[len(right)-1] - left[0])/ len(right))

    # we need a (numrects x numsides x 2) numpy array for the path helper                                                
    # function to build a compound path                                                                                  
    XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T

    # get the Path object                                                                                                
    barpath = path.Path.make_compound_path_from_polys(XY)

    # make a patch out of it                                                                                             
    patch = patches.PathPatch(barpath, facecolor=colour, edgecolor='black', alpha=1)
    ax.add_patch(patch)
    
##    #set major minor axis
##    minorLocatory0 = MultipleLocator(5)
##    ax[0].yaxis.set_minor_locator(minorLocatory0)
    
    #draw legend
    legend_patch = [patches.Patch(color=colour, label = "{:s}".format(r'$e^{-}e^{+}$'r'$\rightarrow$ Z'r'/$\gamma$'r'$\rightarrow$'r'$\mu^{-}\mu^{-}$'))]
    plt.legend(handles=legend_patch, bbox_to_anchor=(0.97, 0.99))
    
    # update the view limits                                                                                             
    ax.set_xlim(array_range[0], array_range[1])
    ax.set_ylim(bottom.min(), top.max()+ 0.08*top.max())
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
##    plt.savefig(Save_path + name + '.pdf')
##    plt.show()

def residual(array, wghts, name, xlabel, ylabel, array_range = 0):
    
    fig, ax = plt.subplots(2,1, sharex=True, tight_layout=True, gridspec_kw = {'height_ratios':[2, 1]})
##    fig, ax = plt.subplots(2,1, sharex=True, tight_layout=True)
    if array_range is 0:
        array_range = (np.min(array), np.max(array))
        
    n, bins = np.histogram(array, 50, weights=wghts)
    print('total CS(pb) = %s' % sum(n))
          
    # get the corners of the rectangles for the histogram                                                                
    left = np.array(bins[:-1])
    right = np.array(bins[1:])
    bottom = np.zeros(len(left))
    width = ((right[len(right)-1] - left[0])/ len(right))
    top = (bottom + n) / width
    A_top = []
    #Calculate analytical
    for i in range(len(left)):
        x1 = right[i]
        x0 = left[i]
        A_top.append(Analytic(x1, x0))    
    # function to build a compound path                                                                                  
    XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T

    # get the Path object                                                                                                
    barpath = path.Path.make_compound_path_from_polys(XY)

    # make a patch out of it                                                                                             
    patch = patches.PathPatch(barpath, facecolor= colour, edgecolor='black', alpha=1)
    ax[0].add_patch(patch)
##    ones = np.ones(len(bins))

    #draw legend
    legend_patch = [patches.Patch(color=colour, label = "{:s}".format(r'$e^{-}e^{+}$'r'$\rightarrow \ Z$'r'/$\gamma$'r'$\rightarrow$'r'$\mu^{-}\mu^{-}$'))]
    ax[0].legend(handles=legend_patch, ncol=1, frameon=False, fontsize= 16)

     #set font dict
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 18,
            }
    
    '''==================================================='''
    '''==================================================='''
    ax[0].set_ylabel(r'$d\sigma / dcos\theta \ [pb]$', fontsize=18)
    #update the view limits
    ax[0].set_xlim(-1, 1)
    ax[0].set_ylim(0, 900)
    #set minor axis
    minorLocator = AutoMinorLocator(5)
    ax[0].xaxis.set_minor_locator(minorLocator)
    minorLocator = AutoMinorLocator(5)
    ax[0].yaxis.set_minor_locator(minorLocator)
    ax[0].tick_params(axis='both', which='major', labelsize=14)
    ax[0].xaxis.labelpad = 0
    ax[0].yaxis.labelpad = 15   
    ax[0].text(-0.70, 780, r'$E_{cm} = 90\ GeV$' , fontdict=font) #r'$math.sqrt(S)=1.96\ TeV$'

    
    '''==================================================='''
    '''==================================================='''
    ax[1].set_xlabel(r'$cos \ \theta$', fontsize=18)
    ax[1].locator_params(tight=True, nbins=6)
##    # update the view limits
##    ax[1].set_xlim(0, 2.25)
##    ax[1].set_ylim(0, 5.2)
    #set minor axis
    minorLocator = AutoMinorLocator(5)
    ax[1].xaxis.set_minor_locator(minorLocator)
    minorLocator = AutoMinorLocator(5)
    ax[1].yaxis.set_minor_locator(minorLocator)
    ax[1].tick_params(axis='both', which='major', labelsize=14)
    ax[1].xaxis.labelpad = 0
    ax[1].yaxis.labelpad = 0    

    
    top1 = top * width    

    ax[1].get_yaxis().get_major_formatter().set_useOffset(False)


    #error
    error = [(x - y)*100/y for x, y in zip(top1, A_top)]
    ax[1].plot(bins[:-1], error, color='red',linewidth=1.8) # plot error
    ax[1].set_ylim(-0.05, 0.05) #error
    ax[1].set_ylabel(r'$Error \ (\%)$', fontsize=18)

    #ratio
##    ratio = [x / y for x, y in zip(top1, A_top)]
##    ax[1].plot(bins[:-1], ratio, color='red',linewidth=1.8) # plot ratio    
##    ax[1].set_ylim(np.min(ratio), np.max(ratio)) #ratio
##    ax[1].set_ylabel('ratio')
##    ax[1].set_ylim(1.0015, 0.9985)

    
    # update the view limits
    bins = [x + width/2 for x in bins]

    
##    ax[0].plot(bins[0:-1], A_top, color='red',linewidth=1.8) # plot analytic
##    ax[0].set_xlim(left[0], right[-1])
##    ax[0].set_ylim(bottom.min(), top.max()+ 0.08*top.max())
##    ax[0].set_ylabel(ylabel)
    
##    plt.savefig(Save_path + name + '.pdf')
    plt.tight_layout()
    plt.savefig(Save_path + name + '.pdf')    
##    plt.show()

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

'----=========================================================================----'
'----==========================      LOAD FILES     ==========================----'
'----=========================================================================----'

cwd = os.getcwd()
np_name = glob.glob('*.npy')[0]
print np_name
load = np.load(np_name)
 #ee_array = [PScosth, PSpT, PSrapidity, wghts]
PScosth = load[0]
PSpT = load[1]
PSrapidity = load[2]
wghts = load[3]
array_names = 'PScosth, PSpT, PSrapidity, wghts'
'----=========================================================================----'
'----==========================      SAVE HISTO     ==========================----'
'----=========================================================================----'

 #makes directories
Save_path = cwd + '/Histo/'
if not os.path.exists(os.path.dirname(Save_path)):
    try:
        os.makedirs(os.path.dirname(Save_path))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise        
## #save text
##with open(os.path.join(Save_path, 'details.txt'), 'w') as myfile:
##    myfile.write('\n\n || Array names = PScosth, PSpT, PSrapidity, wghts ||')

'----=========================================================================----'
'----==========================    HISTOGRAMMING    ==========================----'
'----=========================================================================----'
colour = 'gray'

##print 'Generating cos distribution for ee'
## #generate_histo(PScosth, 'PScosth_Eve', r'cos$\theta$') 
##generate_normalized_histo(PScosth, wghts, 'PScosth', r'cos$\theta$', 'd'r'$\sigma$''/d'r'cos$\theta$')
##
##print 'Generating transverse momentum distribution for ee'
## #generate_histo(PSpT,'PSpT_Eve', r'$P_{T}$') 
##generate_normalized_histo(PSpT, wghts, 'PSpT', r'$P_{T}$', 'd'r'$\sigma$''/d'r'$P_{T}$')

print 'Generating rapidity distribution for ee'
##generate_histo(PSrapidity,'PSrapidity_Eve', r'$y_{e^{-}}$')
##generate_normalized_histo(PSrapidity, wghts, 'PSrapidity', r'$y_{\mu^{-}}$', 'd'r'$\sigma$''/d'r'$y_{\mu^{-}}$', (-5, 5))

print 'Generating residual for costheta'
residual(PScosth, wghts, 'Res_costh', r'cos$\theta$', 'd'r'$\sigma$''/d'r'cos$\theta$ [pb]')
##
##plt.show()
