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
from matplotlib.ticker import AutoMinorLocator

'----=========================================================================----'
'----==========================      LOAD FILES     ==========================----'
'----=========================================================================----'

cwd = os.getcwd()
base = os.path.splitext(os.path.basename(__file__))[0]
aa = []
for np_name in glob.glob('*.npy'):
    array_name = os.path.splitext(np_name)[0]
    load = np.load(np_name)
    aa.append(load)
    print ('%s array loaded!!' % array_name)

soft_CS = np.asarray(aa[0])
tot_hsq_CS = np.asarray(aa[1])
tot_hs_CS = np.asarray(aa[2])
hard_CS = np.asarray(aa[3])
E_cut = np.asarray(aa[4])

 #set font dict
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 18,
        }

fig, ax = plt.subplots()
plt.semilogx(E_cut, soft_CS, label = 'Soft', color='blue', linewidth=2)
plt.semilogx(E_cut, hard_CS, label = 'Hard', color='green', linewidth=2)
plt.semilogx(E_cut, tot_hsq_CS, linestyle='--', label = 'Sum', color='red', linewidth=2)
plt.legend(frameon=False, fontsize=18)
plt.xlabel(r'$E_{cut}$', fontsize=20)
plt.ylabel(r'$\sigma^{tot}_{A}$' r'$(E_{cut})$', fontsize=20)

#change tick label size
ax.tick_params(axis='both', which='major', labelsize=14)
ax.tick_params(axis='both', which='minor', labelsize=8)

ax.xaxis.labelpad = 2
ax.yaxis.labelpad = 0  

plt.text(9.**-3, 1.65, r'$\sqrt{S} = 1.96\ TeV$' , fontdict=font) #r'$math.sqrt(S)=1.96\ TeV$'

axes = plt.gca()
axes.set_xlim([10**-3,E_cut[-1]+10])

plt.savefig('23body.pdf')

plt.tight_layout()
##plt.show()
