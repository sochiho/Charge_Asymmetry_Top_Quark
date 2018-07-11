'''Function for creating and binning Observables'''
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

#comment out the following if not using matplotlib and numpy`
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
import cPickle as pickle
import resource
import gc
from scipy.interpolate import spline
from matplotlib.ticker import AutoMinorLocator

PDF1 = 'cteq66'
PDF2 = 'MSTW2008CPdeutnlo68cl'
PDF3 = 'NNPDF21_lostar_as_0130_100'
PDF_set = [PDF1, PDF2, PDF3]

#Global varibales
ECM = 1960
s = ECM**2
print "hadron com energy:", ECM, "GeV"

#Define variables in Kuhn paper
M_top = 173.1 # Mass of top quark GeV
alpha_QCD = 0.109 # alpha QCD
E_cut = 0.01 * M_top
deltath = 2
tau0 = 4*M_top**2/s

def mem():
    print('Memory usage         : % 8.3f GB' % (
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/10**6))

'----=========================================================================----'
'----=======================     HISTOGRAM BINNING    ========================----'
'----=========================================================================----'

def generate_histo(array, name, xlabel, array_range = 0):
  
    fig, ax = plt.subplots()
    
    #axis limits 
    top_lim= []
    bottom_lim = []
    fcolour_list = ['grey', 'lightgray', 'dimgrey']
    alpha_list = [0.7, 1, 1]
    labelit_list = label_name
    hatch_list = [None, None, None]
    zorder_list = [0, 1, 2]
    #ensure all arrays use same range
    if array_range is 0:
        array_range = (np.min(array[0]), np.max(array[0]))
        print array_range, 'array range'
        
    for ii in range(len(array)):
        n, bins = np.histogram(array[ii], 50, array_range)

        # get the corners of the rectangles for the histogram
        left = np.array(bins[:-1])
        right = np.array(bins[1:])
        bottom = np.zeros(len(left))
        top = bottom + n

        # we need a (numrects x numsides x 2) numpy array for the path helper
        XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T
        # get the Path object
        barpath = path.Path.make_compound_path_from_polys(XY)
        
        # make a patch out of it
        patch = patches.PathPatch(barpath, facecolor= fcolour_list[ii], edgecolor='black', alpha=alpha_list[ii], hatch=hatch_list[ii], zorder = zorder_list[ii])
        ax.add_patch(patch)

        #determine top, bottom limits
        top_lim.append(top.max())
        bottom_lim.append(bottom.min())
        
    #draw legend
    legend_patch = [patches.Patch(color=fcolour_list[i], label = "{:s}".format(labelit_list[i])) for i in range(len(fcolour_list))]
    plt.legend(handles=legend_patch, ncol=1)    
    #set minor axis
    minorLocator = MultipleLocator((array_range[1]-array_range[0])/20)
    ax.xaxis.set_minor_locator(minorLocator)
    # update the view limits
    ax.set_xlim(array_range[0], array_range[1])
    ax.set_ylim(np.min(bottom_lim), np.max(top_lim) + 0.08 * np.max(top_lim))
    
    plt.xlabel(xlabel)
    plt.ylabel('Number of events')

##    plt.savefig(Save_path + name + '.pdf')
##    plt.show()

def generate_normalized_histo(array, wghts, name, xlabel, ylabel, array_range = 0):
    number_bins = 50
    fig, ax = plt.subplots()

    #axis limits 
    top_lim= []
    bottom_lim = []
    fcolour_list = ['dimgray', 'darkgrey', 'gainsboro'] ###
##    fcolour_list = ['b', 'lightskyblue', 'navy']
##    fcolour_list = ['lawngreen', 'turquoise', 'orange']
    alpha_list = [1, 1, 1]
##    labelit_list = [r'$M_{top}$', r'$M_{top}$/2', r'$M_{top}$'r'$\times$ 2'] ###
    labelit_list = label_name
    hatch_list = [None, None, None]
    zorder_list = [0, 1, 2]
##    zorder_list = [0, 2, 1]
    #ensure all arrays use same range
    if array_range is 0:
        array_range = (np.min(array[0]), np.max(array[0]))
    for ii in range(len(array)):
        n, bins = np.histogram(array[ii], number_bins, array_range, weights=wghts[ii])
        width = (bins[-1]- bins[0])/number_bins
        n = n/width
        '========================================================================='
        #smooth out histogram        
        n0 = n[0]
        nn = n[-1]
        binsn = bins[-1]
        n, bins = np.histogram(array[ii], 8, array_range, weights=wghts[ii])
        width = (bins[-1]- bins[0])/8
        n = n/width        

        n = np.append(n,[nn])
        n = np.append([n0], n)
        h_width = (bins[1]-bins[0])/2
        mid_bin = [x + h_width for x in bins[0:-1]]
        mid_bin = np.append(bins[0], mid_bin)
        mid_bin = np.append(mid_bin, [binsn])
        x = np.linspace(bins[0], bins[-1], 50) #smooth back to 50 bins
        n = spline(mid_bin, n, x[0:-1])
        bins = x
        '========================================================================='
        
        print n.max()
        # get the corners of the rectangles for the histogram                                                                        
        left = np.array(bins[:-1])
        right = np.array(bins[1:])
        bottom = np.zeros(len(left))
        top = (bottom + n)#/ ((right[-1] - left[0])/ len(right))

        # we need a (numrects x numsides x 2) numpy array for the path helper                                                                                               
        XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T

        # get the Path object                                                                                                
        barpath = path.Path.make_compound_path_from_polys(XY)

        # make a patch out of it                                                                                             
        patch = patches.PathPatch(barpath, facecolor= fcolour_list[ii], edgecolor='black', alpha=alpha_list[ii], hatch=hatch_list[ii], zorder = zorder_list[ii])
        ax.add_patch(patch)
        
        #determine top, bottom limits
        top_lim.append(top.max())
        bottom_lim.append(bottom.min())

     #set font dict
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 16,
            }
    #draw legend
    legend_patch = [patches.Patch(color=fcolour_list[i], label = "{:s}".format(labelit_list[i])) for i in range(len(labelit_list))]
    plt.legend(handles=legend_patch, ncol=1, frameon=False, fontsize = 16)
    
    #set minor x axis
    minorLocator = AutoMinorLocator(5)
    ax.yaxis.set_minor_locator(minorLocator)
    
    # update the view limits
    ax.set_xlim(array_range[0], array_range[1])
    ax.set_ylim(np.min(bottom_lim), np.max(top_lim) + 0.08 * np.max(top_lim))
    #change tick label size
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='minor', labelsize=8)    
    #changge x,y fontsize
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=20)

    #print Tevatron power
    xp = array_range[0] + (array_range[1] - array_range[0]) / 3
    yp = np.max(top_lim) - np.max(top_lim)/20
    plt.text(xp, yp, r'$\sqrt{S}=1.96\ TeV$' , fontdict=font) #r'$math.sqrt(S)=1.96\ TeV$'
    
    plt.savefig(Save_path + name + '.pdf')
    plt.tight_layout()
    plt.show()

def generate_asymm_plot(array, array_LO, wghts, wghts_LO, name, xlabel, ylabel, array_range = 0):

    fig, ax = plt.subplots()

    #axis limits 
    top_lim= []
##    fcolour_list = ['grey', 'lightgray', 'dimgrey'] ###
##    fcolour_list = ['b', 'lightskyblue', 'navy']
    fcolour_list = ['lawngreen', 'turquoise', 'orange']
    alpha_list = [1, 1, 1]
##    labelit_list = [r'$M_{top}$', r'$M_{top}$/2', r'$M_{top}$'r'$\times$ 2'] ###
    labelit_list = label_name
    hatch_list = [None, None, None]
##    zorder_list = [1, 0, 2]
    zorder_list = [0, 1, 2]
    #ensure all arrays use same range
    if array_range is 0:
        array_range = (np.min(array[0]), np.max(array[0]))

    plots = []

    for ii in range(len(array)):
        number_bins = 6
        n, bins = np.histogram(array[ii], number_bins, array_range, weights=wghts[ii])
        n_LO, bins_LO = np.histogram(array_LO[ii], number_bins, array_range, weights=wghts_LO[ii])
        
        width = (bins[-1]- bins[0])/number_bins
        n = n/width
        n_LO = n_LO/width

        '========================================================================='
        top = [0]
        for i in range(len(n)):
            if n_LO[i] == 0:
                t = 0
            else:
                t = n[i]/n_LO[i]
            top.append(t)
        plots.append(top)
        
        #determine top, bottom limits
        top_lim.append(np.max(top))
        
        '========================================================================='
     #set font dict
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 16,
            }
    #smooth out curve
    xnew = np.linspace(bins[0], bins[-1],300) #300 represents number of points to make between T.min and T.max
    plots[0] = spline(bins, plots[0], xnew)
    plots[1] = spline(bins, plots[1], xnew)
    plots[2] = spline(bins, plots[2], xnew)         
    
    bins = xnew
    '========================================================================='

    #draw legend
##    legend_patch = [patches.Patch(color=fcolour_list[i], label = "{:s}".format(labelit_list[i])) for i in range(len(labelit_list))]
##    plt.legend(handles=legend_patch, bbox_to_anchor=(0.285, 0.99), frameon=False, fontsize = 16)

    #Convert plots to arrays
    plots[0] = np.asarray(plots[0])
    plots[1] = np.asarray(plots[1])
    plots[2] = np.asarray(plots[2])

    ax.plot(bins, plots[0], color='black', label= labelit_list[0], linewidth=2.5, alpha=1, linestyle= '-')
    ax.plot(bins, plots[1], color='black', label= labelit_list[1], linewidth=1.6, alpha=1, linestyle= '--')
    ax.plot(bins, plots[2], color='black', label= labelit_list[2], linewidth=0.9, alpha=1, linestyle= '-')

    #draw legend
    plt.legend(bbox_to_anchor=(0.3, 0.97), frameon= False, fontsize = 16)
    
    ax.fill_between(bins, plots[0], plots[2], where=plots[0] >= plots[2], facecolor='grey', interpolate=True)    
##    ax.fill_between(bins, plots[1], plots[2], where=plots[1] >= plots[2], facecolor='grey', interpolate=True)    

    #set minor x axis
    minorLocator = AutoMinorLocator(5)
    ax.xaxis.set_minor_locator(minorLocator)
    minorLocator = AutoMinorLocator(5)
    ax.yaxis.set_minor_locator(minorLocator)
                   
    # update the view limits
    ax.set_xlim(array_range[0], array_range[1])
    ax.set_ylim(0, np.max(top_lim) + 0.08 * np.max(top_lim))
    #change tick label size
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='minor', labelsize=8)    
    #changge x,y fontsize
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=20)

    #print Tevatron power
    xp = array_range[0] + (array_range[1] - array_range[0]) / 3
    yp = np.max(top_lim) - np.max(top_lim)/20
    plt.text(xp, yp, r'$\sqrt{S}=1.96\ TeV$' , fontdict=font) #r'$math.sqrt(S)=1.96\ TeV$'

    plt.savefig(Save_path + name + '.pdf')
    plt.tight_layout()    
    plt.show()

'----=========================================================================----'
'----=========================      DICTIONARIES     =========================----'
'----=========================================================================----'

#Create observable dictionaries

Dict = {}

## #key finding
##for key in Dict:
##    print key, 'PDF key\n'
##    for kk in Dict[key]:
##        print kk, 'Observables'
##
##print Dict['MS']['y3_LO']

'----=========================================================================----'
'----=========================================================================----'
'----=========================      OBSERVABLES     ==========================----'
'----=========================================================================----'
'----=========================================================================----'

cwd = os.getcwd()
base = os.path.splitext(os.path.basename(__file__))[0]

'----=========================================================================----'

#Impost cut either 'ON' or 'OFF'
CUT = 'ON_y'

'----=========================================================================----'

 #Frame = 'LAB' or 'tt'
if 'LAB' in cwd:
    FRAME = 'LAB'
    print '%s frame' % FRAME
elif 'tt' in cwd:
    FRAME = 'tt'
    print '%s frame' % FRAME

'----=========================================================================----'
'----==========================      LOAD FILES     ==========================----'
'----=========================================================================----'
print time.strftime("%H:%M")

for pickle_name in glob.glob('*.pickle'):
    mem()
    '''Load pickle file and produce observables'''
    dict_name = os.path.splitext(pickle_name)[0]
    print ('%s dictionary loading!!\n' % dict_name)
    ss = dict_name.split('_')
    mu_f = ss[1]+'_'+ss[2]
    int_region = ss[3]
    PDF = ss[0]
    
    #create PDF and mu_f dictionary in Dict
    if Dict.get(PDF) == None:
        Dict[PDF] = {}
    d_PDF = Dict[PDF] # PDF dictionary path for Dict
    if d_PDF.get(mu_f) == None:
        d_PDF[mu_f] = {}
    d_mu = d_PDF[mu_f]
    
    pickle_in = open(pickle_name, 'rb')
    data = pickle.load(pickle_in)

    '----=========================================================================----'
    '----========================     qq CALCULATION    ==========================----'
    '----=========================================================================----'
    counter = 1
    if int_region == 'qqbar':
        qq_bar = data
        print '----=========================================================================---- \n'
        print ('Generating qq_bar observables, PDF = %s, mu = %s \n' % (PDF, mu_f))
        #Random numbers from VEGAS
        R0 = qq_bar.get('x0')
        R1 = qq_bar.get('x1')
        R2 = qq_bar.get('x2')
        wght = qq_bar.get('wghts')

        jj=0
        Neve = len(R0)
        for ii in range(Neve):
            sys.stdout.write("progress: %d%%   \r" % (float(jj)*100./(Neve)) )
            sys.stdout.flush()            
            x1 = tau0 + (1 -tau0) * R0[ii]
            x2 = tau0/x1 +(1-tau0/x1)*R1[ii]
            hats = x1*x2*s
            xm = M_top**2/hats
            Q = math.sqrt(hats)
            costh = -1 + deltath*R2[ii]
            E3 = 0.5
            k3 = math.sqrt(E3**2 -xm)
            y3 = np.log((E3 + k3*costh)/(E3 - k3*costh))/2
            y4 = np.log((E3 - k3*costh)/(E3 + k3*costh))/2
            ydelta = y3-y4
            if FRAME == 'LAB':
                y3 = y3 + np.log(x1/x2)/2
            elif FRAME == 'tt':
                y3 = y3

            #append to LO list
            if FRAME == 'LAB':
                add = {'y3_LO':y3, 'wghts_LO':wght[ii]}                
            
            if FRAME == 'tt':
                if CUT == 'ON':
                    if Q < 450:
                        add = {'Qs_LO':Q, 'swghts_LO':wght[ii], 'swghts2_LO': wght[ii]*1000}
                    else:
                        add = {'Ql_LO':Q, 'lwghts_LO':wght[ii], 'lwghts2_LO': wght[ii]*1000}
                elif CUT == 'ON_y':
                    if abs(ydelta) < 1:
                        add = {'Qs_LO':Q, 'ydeltas_LO': ydelta, 'swghts_LO':wght[ii], 'swghts2_LO': wght[ii]*1000}
                    else:
                        add = {'Ql_LO':Q, 'ydeltal_LO': ydelta, 'lwghts_LO':wght[ii], 'lwghts2_LO': wght[ii]*1000}
                        
                else:
                    add = {'y3_LO':y3, 'ydelta_LO':ydelta, 'Q_LO':Q, 'wghts_LO':wght[ii], 'wghts2_LO': wght[ii]*1000}
            
            for key, val in add.items():
                if key in d_mu:
                    if counter in range(1,5):
                        print 'appending to %s key in %s' % (key, PDF)
                    d_mu[key].append(val)
                else:
                    print 'new mu %s key created in %s!!' % (key, PDF)
                    d_mu[key] = []
                    d_mu[key].append(val)
                counter +=1
            jj +=1
            
        print '\n=================================='   
        print ' qq_bar CS(pb) = %s' % sum(wght)
        print '==================================\n'   

    '----=========================================================================----'
    '----========================     gg CALCULATION    ==========================----'
    '----=========================================================================----'

    counter = 1
    if int_region == 'gg':
        gg = data 
        print '----=========================================================================---- \n'
        print ('Generating gg observables, PDF = %s, mu = %s \n' % (PDF, mu_f))
        #Random numbers from VEGAS
        R0 = gg.get('x0')
        R1 = gg.get('x1')
        R2 = gg.get('x2')
        wght = gg.get('wghts')
        #Observable arrays
        y3 = []
        y3lab = []
        ydelta = []
        Q = []
        
        jj=0
        Neve = len(R0)
        for ii in range(Neve):
            sys.stdout.write("progress: %d%%   \r" % (float(jj)*100./(Neve)) )
            sys.stdout.flush()            
            x1 = tau0 + (1 -tau0) * R0[ii]
            x2 = tau0/x1 +(1-tau0/x1)*R1[ii]
            hats = x1*x2*s
            xm = M_top**2/hats
            Q = math.sqrt(hats)
            costh = -1 + deltath*R2[ii]
            E3 = 0.5
            k3 = math.sqrt(E3**2 -xm)
            y3 = np.log((E3 + k3*costh)/(E3 - k3*costh))/2
            y4 = np.log((E3 - k3*costh)/(E3 + k3*costh))/2
            ydelta = y3-y4
            if FRAME == 'LAB':
                y3 = y3 + np.log(x1/x2)/2
            elif FRAME == 'tt':
                y3 = y3
                                
            #append to LO list
            if FRAME == 'LAB':
                add = {'y3_LO':y3, 'wghts_LO':wght[ii]}                
            
            if FRAME == 'tt':
                if CUT == 'ON':
                    if Q < 450:
                        add = {'Qs_LO':Q, 'swghts_LO':wght[ii], 'swghts2_LO': wght[ii]*1000}
                    else:
                        add = {'Ql_LO':Q, 'lwghts_LO':wght[ii], 'lwghts2_LO': wght[ii]*1000}
                        
                elif CUT == 'ON_y':
                    if abs(ydelta) < 1:
                        add = {'Qs_LO':Q, 'ydeltas_LO': ydelta, 'swghts_LO':wght[ii], 'swghts2_LO': wght[ii]*1000}
                    else:
                        add = {'Ql_LO':Q, 'ydeltal_LO': ydelta, 'lwghts_LO':wght[ii], 'lwghts2_LO': wght[ii]*1000}
                        
                else:
                    add = {'y3_LO':y3, 'ydelta_LO':ydelta, 'Q_LO':Q, 'wghts_LO':wght[ii], 'wghts2_LO': wght[ii]*1000}
            
            for key, val in add.items():
                if key in d_mu:
                    if counter in range(1,5):
                        print 'appending to %s key in %s' % (key, PDF)
                    d_mu[key].append(val)
                else:
                    print 'new mu %s key created in %s!!' % (key, PDF)
                    d_mu[key] = []
                    d_mu[key].append(val)
                counter +=1
            jj +=1
            
        print '\n=================================='                       
        print '  gg CS(pb) = %s' % sum(wght)
        print '==================================\n'   
        

    '----=========================================================================----'
    '----========================       SOFT REGION     ==========================----'
    '----=========================================================================----'

    counter = 1
    if int_region == 'soft':
        soft = data 
        print '----=========================================================================---- \n'
        print ('Generating Soft observables, PDF = %s, mu = %s \n' % (PDF, mu_f))
        #Random numbers from VEGAS
        R0 = soft.get('x0')
        R1 = soft.get('x1')
        R2 = soft.get('x2')
        wght = soft.get('wghts')

        jj=0
        Neve = len(R0)
        for ii in range(Neve):
            sys.stdout.write("progress: %d%%   \r" % (float(jj)*100./(Neve)) )
            sys.stdout.flush()            
            x1 = tau0 + (1 -tau0) * R0[ii]
            x2 = tau0/x1 +(1-tau0/x1)*R1[ii]
            hats = x1*x2*s
            xm = M_top**2/hats
            Q = math.sqrt(hats)
            costh = -1 + deltath*R2[ii]
            E3 = 0.5
            k3 = math.sqrt(E3**2 -xm)
            y3 = np.log((E3 + k3*costh)/(E3 - k3*costh))/2
            y4 = np.log((E3 - k3*costh)/(E3 + k3*costh))/2
            ydelta = y3-y4
            if FRAME == 'LAB':
                y3 = y3 + np.log(x1/x2)/2
            elif FRAME == 'tt':
                y3 = y3

            #append to LO list
            if FRAME == 'LAB':
                add = {'y3_AS':y3, 'wghts_AS':wght[ii]}                
            
            if FRAME == 'tt':
                if CUT == 'ON':
                    if Q < 450:
                        add = {'Qs_AS':Q, 'swghts_AS':wght[ii], 'swghts2_AS': wght[ii]*1000}
                    else:
                        add = {'Ql_AS':Q, 'lwghts_AS':wght[ii], 'lwghts2_AS': wght[ii]*1000}

                elif CUT == 'ON_y':
                    if abs(ydelta) < 1:
                        add = {'Qs_AS':Q, 'ydeltas_AS': ydelta, 'swghts_AS':wght[ii], 'swghts2_AS': wght[ii]*1000}
                    else:
                        add = {'Ql_AS':Q, 'ydeltal_AS': ydelta, 'lwghts_AS':wght[ii], 'lwghts2_AS': wght[ii]*1000}
                        
                else:
                    add = {'y3_AS':y3, 'ydelta_AS':ydelta, 'Q_AS':Q, 'wghts_AS':wght[ii], 'wghts2_AS': wght[ii]*1000}

            for key, val in add.items():
                if key in d_mu:
                    if counter in range(1,5):
                        print 'appending to %s key in %s' % (key, PDF)
                    d_mu[key].append(val)
                else:
                    print 'new mu %s key created in %s!!' % (key, PDF)
                    d_mu[key] = []
                    d_mu[key].append(val)
                counter +=1
            jj +=1
            
        print '\n=================================='   
        print '  soft CS(pb) = %s' % sum(wght)
        print '==================================\n'   
    
    '----=========================================================================----'
    '----========================       HARD REGION     ==========================----'
    '----=========================================================================----'

    counter = 1
    if int_region == 'hard':
        hard = data 
        print '----=========================================================================---- \n'
        print ('Generating Hard observables, PDF = %s, mu = %s \n' % (PDF, mu_f))
        #Random numbers from VEGAS
        R0 = hard.get('x0')
        R1 = hard.get('x1')
        R2 = hard.get('x2')
        R3 = hard.get('x3')
        R4 = hard.get('x4')
        R5 = hard.get('x5')
        wght = hard.get('wghts')

        jj=0
        Neve = len(R0)
        for ii in range(Neve):
            sys.stdout.write("progress: %d%%   \r" % (float(jj)*100./(Neve)) )
            sys.stdout.flush()            
            x1 = tau0 + (1 -tau0) * R0[ii]
            x2 = tau0/x1 +(1-tau0/x1)*R1[ii]
            hats = x1*x2*s
            xm = M_top**2/hats
            Q = math.sqrt(hats)
            b = math.sqrt(1-4*xm)
            w = E_cut/math.sqrt(hats)   

            costh = -1 + deltath*R2[ii]
            sinth = math.sqrt(1-costh**2)
            phi = 2*math.pi*R3[ii]
            cosphi = np.cos(phi)
            sinphi = np.sin(phi)

            y45_min = w*(1-b)
            y45_max = 1-2*math.sqrt(xm)
            if y45_max-y45_min < 0:
                continue
            y45 = y45_min + (y45_max-y45_min)*R4[ii]

            E3 = (1-y45)/2
            k3 = math.sqrt(E3**2-xm)
            
            y35_min = max((E3-k3)*y45/(1-E3+k3), 2*w-y45)
            y35_max = (E3+k3)*y45/(1-E3-k3)
            y35 = y35_min + (y35_max-y35_min)*R5[ii]

            E5 = (y35+y45)/2
            if E5< w:
                continue
            cosa = (y35/(2*E5) -E3) /k3
            sina = math.sqrt(1-cosa**2)
            
            E4 = 1-E3-E5
            k4z = -E5*sina*sinth*cosphi + (-k3+E5*cosa)*costh
            y3 = np.log((E3+k3*costh)/(E3-k3*costh)) / 2
            y4 = np.log((E4+k4z)/(E4-k4z)) / 2
            ydelta = y3- y4
            if FRAME == 'LAB': 
                y3 = y3 + np.log(x1/x2)/2
            elif FRAME == 'tt':
                y3 = y3
                
            #append to LO list
            if FRAME == 'LAB':
                add = {'y3_AS':y3, 'wghts_AS':wght[ii]}                
            
            if FRAME == 'tt':
                if CUT == 'ON':
                    if Q < 450:
                        add = {'Qs_AS':Q, 'swghts_AS':wght[ii], 'swghts2_AS': wght[ii]*1000}
                    else:
                        add = {'Ql_AS':Q, 'lwghts_AS':wght[ii], 'lwghts2_AS': wght[ii]*1000}

                elif CUT == 'ON_y':
                    if abs(ydelta) < 1:
                        add = {'Qs_AS':Q, 'ydeltas_AS': ydelta, 'swghts_AS':wght[ii], 'swghts2_AS': wght[ii]*1000}
                    else:
                        add = {'Ql_AS':Q, 'ydeltal_AS': ydelta, 'lwghts_AS':wght[ii], 'lwghts2_AS': wght[ii]*1000}
                        
                else:
                    add = {'y3_AS':y3, 'ydelta_AS':ydelta, 'Q_AS':Q, 'wghts_AS':wght[ii], 'wghts2_AS': wght[ii]*1000}

            for key, val in add.items():
                if key in d_mu:
                    if counter in range(1,5):
                        print 'appending to %s key in %s' % (key, PDF)
                    d_mu[key].append(val)
                else:
                    print 'new mu %s key created in %s!!' % (key, PDF)
                    d_mu[key] = []
                    d_mu[key].append(val)
                counter +=1
            jj +=1

        print '\n=================================='   
        print '  hard CS(pb) = %s' % sum(wght)
        print '==================================\n'   


        '----=========================================================================----'
    '----========================        QG REGION      ==========================----'
    '----=========================================================================----'

    counter = 1
    if int_region == 'qg':
        qg = data
        print '----=========================================================================---- \n'
        print ('Generating qg observables, PDF = %s, mu = %s \n' % (PDF, mu_f))
        #Random numbers from VEGAS
        R0 = qg.get('x0')
        R1 = qg.get('x1')
        R2 = qg.get('x2')
        R3 = qg.get('x3')
        R4 = qg.get('x4')
        R5 = qg.get('x5')
        wght = qg.get('wghts')

        jj=0
        Neve = len(R0)
        for ii in range(Neve):
            sys.stdout.write("progress: %d%%   \r" % (float(jj)*100./(Neve)) )
            sys.stdout.flush()            
            x1 = tau0 + (1 -tau0) * R0[ii]
            x2 = tau0/x1 +(1-tau0/x1)*R1[ii]
            hats = x1*x2*s
            xm = M_top**2/hats
            Q = math.sqrt(hats)
            b = math.sqrt(1-4*xm)
            w = E_cut/math.sqrt(hats)   

            costh = -1 + deltath*R2[ii]
            sinth = math.sqrt(1-costh**2)
            phi = 2*math.pi*R3[ii]
            cosphi = np.cos(phi)
            sinphi = np.sin(phi)

            y45_min = w*(1-b)
            y45_max = 1-2*math.sqrt(xm)
            if y45_max-y45_min < 0:
                continue
            y45 = y45_min + (y45_max-y45_min)*R4[ii]

            E3 = (1-y45)/2
            k3 = math.sqrt(E3**2-xm)
            
            y35_min = max((E3-k3)*y45/(1-E3+k3), 2*w-y45)
            y35_max = (E3+k3)*y45/(1-E3-k3)
            y35 = y35_min + (y35_max-y35_min)*R5[ii]

            E5 = (y35+y45)/2
            if E5< w:
                continue
            cosa = (y35/(2*E5) -E3) /k3
            sina = math.sqrt(1-cosa**2)
            
            E4 = 1-E3-E5
            k4z = -E5*sina*sinth*cosphi + (-k3+E5*cosa)*costh
            #qg
            y3qg = np.log((E3+k3*costh)/(E3-k3*costh)) / 2
            y4qg = np.log((E4+k4z)/(E4-k4z)) / 2
            yqg_delta = y3qg - y4qg
            
            #switch between lab and tt
            if FRAME == 'LAB':
                y3qg = y3qg + np.log(x1/x2)/2
            elif FRAME == 'tt':
                y3qg = y3qg
                
            #append to AS list
            
            if y3qg > 0:
                
            #append to LO list
                if FRAME == 'LAB':
                    add = {'y3_AS':y3qg, 'wghts_AS':wght[ii]}                
                
                if FRAME == 'tt':
                    if CUT == 'ON':
                        if Q < 450:
                            add = {'Qs_AS':Q, 'swghts_AS':wght[ii], 'swghts2_AS': wght[ii]*1000}
                        else:
                            add = {'Ql_AS':Q, 'lwghts_AS':wght[ii], 'lwghts2_AS': wght[ii]*1000}

                    elif CUT == 'ON_y':
                        if abs(yqg_delta) < 1:
                            add = {'Qs_AS':Q, 'ydeltas_AS': yqg_delta, 'swghts_AS':wght[ii], 'swghts2_AS': wght[ii]*1000}
                        else:
                            add = {'Ql_AS':Q, 'ydeltal_AS': yqg_delta, 'lwghts_AS':wght[ii], 'lwghts2_AS': wght[ii]*1000}
                        
                    else:
                        add = {'y3_AS':y3qg, 'ydelta_AS':yqg_delta, 'Q_AS':Q, 'wghts_AS':wght[ii], 'wghts2_AS': wght[ii]*1000}

                for key, val in add.items():
                    if key in d_mu:
                        if counter in range(1,5):
                            print 'appending to %s key in %s' % (key, PDF)
                        d_mu[key].append(val)
                    else:
                        print 'new mu %s key created in %s!!' % (key, PDF)
                        d_mu[key] = []
                        d_mu[key].append(val)
                    counter +=1
                
            #gq
            y3gq = np.log((E3 - k3*costh)/(E3 + k3*costh)) / 2
            y4gq = np.log((E4 - k4z)/(E4 + k4z)) / 2
            ygq_delta = y3gq-y4gq
            
            #switch between lab and tt
            if FRAME == 'LAB':
                y3gq = y3gq + np.log(x1/x2)/2
            elif FRAME == 'tt':
                y3gq = y3gq
                
            if y3gq > 0:

            #append to LO list
                if FRAME == 'LAB':
                    add = {'y3_AS':y3gq, 'wghts_AS':wght[ii]}                
                
                if FRAME == 'tt':
                    if CUT == 'ON':
                        if Q < 450:
                            add = {'Qs_AS':Q, 'swghts_AS':wght[ii], 'swghts2_AS': wght[ii]*1000}
                        else:
                            add = {'Ql_AS':Q, 'lwghts_AS':wght[ii], 'lwghts2_AS': wght[ii]*1000}

                    elif CUT == 'ON_y':
                        if abs(ygq_delta) < 1:
                            add = {'Qs_AS':Q, 'ydeltas_AS': ygq_delta, 'swghts_AS':wght[ii], 'swghts2_AS': wght[ii]*1000}
                        else:
                            add = {'Ql_AS':Q, 'ydeltal_AS': ygq_delta, 'lwghts_AS':wght[ii], 'lwghts2_AS': wght[ii]*1000}
                            
                    else:
                        add = {'y3_AS':y3gq, 'ydelta_AS':ygq_delta, 'Q_AS':Q, 'wghts_AS':wght[ii], 'wghts2_AS': wght[ii]*1000}                        

                for key, val in add.items():
                    if key in d_mu:
                        if counter in range(1,5):
                            print 'appending to %s key in %s' % (key, PDF)
                        d_mu[key].append(val)
                    else:
                        print 'new mu %s key created in %s!!' % (key, PDF)
                        d_mu[key] = []
                        d_mu[key].append(val)
                    counter +=1
            jj+=1

        print '\n=================================='   
        print 'quark-gluon CS(pb) = %s' % sum(wght)
        print '==================================\n'

        pickle_in.close()
        gc.collect()
        mem()

 #key finding
for key in Dict:
    for key1 in Dict[key]:
        print '\nPDF: %s || mu_f: %s' % (key, key1)
        print '======================'   
        for key2 in Dict[key][key1]:
            print key2, 'Observables'
        print '======================\n'   

'----=========================================================================----'   
'----=========================================================================----'
'----===========================    FB ASYMMETRY    ==========================----'
'----=========================================================================----'
'----=========================================================================----'

#PDF sets
CTEQ = 'cteq66'
MSTW_NLO = 'MSTW2008nlo90cl'
NNPDF = 'NNPDF21_lo_as_0119_100'
##NNPDF_STAR = 'NNPDF21_lostar_as_0119_100'
PDF_set = [CTEQ, MSTW_NLO, NNPDF]

print 'Calculating asymmetries\n'

if FRAME == 'tt' and CUT == 'ON':
    with open('FB_CUT.txt', 'w') as myfile:
        myfile.write('\n FORWARD-BACKWARD ASYMMETRY (%s frame, CUT = %s)' % (FRAME, CUT))                        
    for PDF_key in Dict:
        PDF = [x for x in PDF_set if PDF_key in x][0]
        for mu_key in Dict[PDF_key]:
            print '\nPDF: %s || mu_f: %s' % (PDF_key, mu_key)
            ref = Dict[PDF_key][mu_key]
            # M < 450 GeV
            sNEVAL = len(ref['swghts_LO'])
            sLO_CS = np.sum(ref['swghts_LO'])
            sA_CS = np.sum(ref['swghts_AS'])
            sFB_Asymm = sA_CS / sLO_CS *100
            
            print 'A_CS: %s (M < 450 GeV)' % sA_CS
            print 'LO_CS: %s (M < 450 GeV)' % sLO_CS
            print 'FB: %s %% (M < 450 GeV)' % sFB_Asymm
            print 'NEVAL: %s (M < 450 GeV)' % sNEVAL
            print 'PDF set: %s' % PDF
            
            # M > 450 GeV
            lNEVAL = len(ref['lwghts_LO'])
            lLO_CS = np.sum(ref['lwghts_LO'])
            lA_CS = np.sum(ref['lwghts_AS'])
            lFB_Asymm = lA_CS / lLO_CS *100
            
            print 'A_CS: %s (M > 450 GeV)' % lA_CS
            print 'LO_CS: %s (M > 450 GeV)' % lLO_CS
            print 'FB: %s %% (M > 450 GeV)' % lFB_Asymm
            print 'NEVAL: %s (M > 450 GeV)' % lNEVAL
            print 'PDF set: %s' % PDF
            
            with open('FB_CUT.txt', 'a') as myfile:
                myfile.write('\n\n============================\n')
                myfile.write('\n PDF: %s || mu_f: %s (M < 450 GeV)\n\n' % (PDF_key, mu_key))
                myfile.write('|| A_CS: %s\n|| LO_CS: %s\n|| FB: %s%%\n|| NEVAL: %s\n|| PDF set: %s\n' % (sA_CS, sLO_CS, sFB_Asymm, sNEVAL, PDF))
                myfile.write('\n PDF: %s || mu_f: %s (M > 450 GeV)\n\n' % (PDF_key, mu_key))
                myfile.write('|| A_CS: %s\n|| LO_CS: %s\n|| FB: %s%%\n|| NEVAL: %s\n|| PDF set: %s\n' % (lA_CS, lLO_CS, lFB_Asymm, lNEVAL, PDF))
                myfile.write('============================\n')

if FRAME == 'tt' and CUT == 'ON_y':
    with open('FB_yCUT.txt', 'w') as myfile:
        myfile.write('\n FORWARD-BACKWARD ASYMMETRY (%s frame, CUT = %s)' % (FRAME, CUT))                        
    for PDF_key in Dict:
        PDF = [x for x in PDF_set if PDF_key in x][0]
        for mu_key in Dict[PDF_key]:
            print '\nPDF: %s || mu_f: %s' % (PDF_key, mu_key)
            ref = Dict[PDF_key][mu_key]
            # M < 450 GeV
            sNEVAL = len(ref['swghts_LO'])
            sLO_CS = np.sum(ref['swghts_LO'])
            sA_CS = np.sum(ref['swghts_AS'])
            sFB_Asymm = sA_CS / sLO_CS *100
            
            print 'A_CS: %s (M < 450 GeV)' % sA_CS
            print 'LO_CS: %s (M < 450 GeV)' % sLO_CS
            print 'FB: %s %% (M < 450 GeV)' % sFB_Asymm
            print 'NEVAL: %s (M < 450 GeV)' % sNEVAL
            print 'PDF set: %s' % PDF
            
            # M > 450 GeV
            lNEVAL = len(ref['lwghts_LO'])
            lLO_CS = np.sum(ref['lwghts_LO'])
            lA_CS = np.sum(ref['lwghts_AS'])
            lFB_Asymm = lA_CS / lLO_CS *100
            
            print 'A_CS: %s (M > 450 GeV)' % lA_CS
            print 'LO_CS: %s (M > 450 GeV)' % lLO_CS
            print 'FB: %s %% (M > 450 GeV)' % lFB_Asymm
            print 'NEVAL: %s (M > 450 GeV)' % lNEVAL
            print 'PDF set: %s' % PDF
            
            with open('FB_yCUT.txt', 'a') as myfile:
                myfile.write('\n\n============================\n')
                myfile.write('\n PDF: %s || mu_f: %s (M < 450 GeV)\n\n' % (PDF_key, mu_key))
                myfile.write('|| A_CS: %s\n|| LO_CS: %s\n|| FB: %s%%\n|| NEVAL: %s\n|| PDF set: %s\n' % (sA_CS, sLO_CS, sFB_Asymm, sNEVAL, PDF))
                myfile.write('\n PDF: %s || mu_f: %s (M > 450 GeV)\n\n' % (PDF_key, mu_key))
                myfile.write('|| A_CS: %s\n|| LO_CS: %s\n|| FB: %s%%\n|| NEVAL: %s\n|| PDF set: %s\n' % (lA_CS, lLO_CS, lFB_Asymm, lNEVAL, PDF))
                myfile.write('|| max lower Q_LO: %s\n|| max lower Q_AS: %s\n|| min upper Q_LO: %s\n|| min upper Q_AS: %s\n' % (np.max(ref['Qs_LO']), np.max(ref['Qs_AS']), np.max(ref['Ql_LO']), np.max(ref['Ql_AS'])))
                myfile.write('============================\n')
                
                
else:
    with open('FB.txt', 'w') as myfile:
        myfile.write('\n FORWARD-BACKWARD ASYMMETRY (%s frame, CUT = %s)' % (FRAME, CUT))    
    for PDF_key in Dict:
        PDF = [x for x in PDF_set if PDF_key in x][0]
        for mu_key in Dict[PDF_key]:
            print '\nPDF: %s || mu_f: %s' % (PDF_key, mu_key)
            ref = Dict[PDF_key][mu_key]
            NEVAL = len(ref['wghts_LO'])
            LO_CS = np.sum(ref['wghts_LO'])
            A_CS = np.sum(ref['wghts_AS'])
            FB_Asymm = A_CS / LO_CS *100
            print 'A_CS: %s' % A_CS
            print 'LO_CS: %s' % LO_CS
            print 'FB: %s %%' % FB_Asymm
            print 'NEVAL: %s' % NEVAL
            print 'PDF set: %s' % PDF
            
            with open('FB.txt', 'a') as myfile:
                myfile.write('\n\n============================\n')
                myfile.write('\n PDF: %s || mu_f: %s\n\n' % (PDF_key, mu_key))
                myfile.write('|| A_CS: %s\n|| LO_CS: %s\n|| FB: %s %%\n|| NEVAL: %s\n|| PDF set: %s\n' % (A_CS, LO_CS, FB_Asymm, NEVAL, PDF))
                myfile.write('============================\n')


'----=========================================================================----'
'----==========================      SAVE HISTO     ==========================----'
'----=========================================================================----'
 #make directory path

if len(Dict) == 1:
    s = [key for key in Dict]
    Save_path = cwd+'/'+s[0]+'Histo/'
else:
    Save_path = cwd+'/PDF_Histo/'

 #makes directories
if not os.path.exists(os.path.dirname(Save_path)):
    try:
        os.makedirs(os.path.dirname(Save_path))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

'----=========================================================================----'
'----=========================================================================----'
'----=========================================================================----'

PDF_key = []
for PDF_keys in Dict:
    PDF_key.append(PDF_keys)
    
'----=========================================================================----'
'----=========================================================================----'

mu_key = []
if len(PDF_key) == 1:
    for mu_keys in Dict[PDF_key[0]]:
        mu_key.append(PDF_key[0] +'_' +mu_keys)

    #order arrays
    for i in range(len(PDF_key)):
        order = ['M_half', 'M_topp', 'M_doub']
##        print mu_key, '1st mu_key'
        keys_ordered = []
        label_name = []
        j = 0

        if len(mu_key) != 3:
            keys_ordered = mu_key
            
        else:
            for ii in range(len(order)):
                name = [x for x in mu_key if order[ii] in x][0]
                keys_ordered.append(name)
                if 'M_half' in name:
                    label_name.append(r'$M_{top}/$ $2$')
                if 'M_topp' in name:
                    label_name.append(r'$M_{top}$')
                if 'M_doub' in name:
                    label_name.append(r'$M_{top}$'r'$\times$ $2$')
        
        mu_key = keys_ordered
##        print mu_key, '2nd mu_key'
        
     #produce label names for legend
        
'----=========================================================================----'
'----=========================================================================----' 

if len(PDF_key) != 1:
    i = 0
    for PDF_keys in Dict:
        for mu_keys in Dict[PDF_keys]:
            mu_key.append(PDF_key[i] +'_'+mu_keys)
        i +=1
    l_order = [r'MSTW2008', r'cteq66', r'NNPDF21']
    label_name = []
    
    for i in range(len(PDF_key)):      
        PDF = [x for x in l_order if PDF_key[i] in x][0]
        label_name.append(PDF)
        
'----=========================================================================----'
'----==========================      ARRAY HISTO     =========================----'
'----=========================================================================----'

wghts_LO = []
wghts2_LO = []
y3_LO = []
ydelta_LO = []
Q_LO = []
wghts_AS = []
wghts2_AS = []
y3_AS = []
ydelta_AS = []
Q_AS = []

if CUT != 'ON' or 'ON_y':
    for i in range(len(mu_key)):
        ss = mu_key[i].split('_')
        PDF = ss[0]
        mu = ss[1]+'_'+ss[2]
        Obs_key = Dict[PDF][mu]
        #append Observables
        wghts_AS.append(Obs_key['wghts_AS'])
        wghts_LO.append(Obs_key['wghts_LO'])
        if FRAME == 'LAB':
            y3_AS.append(Obs_key['y3_AS'])
            y3_LO.append(Obs_key['y3_LO'])
        if FRAME == 'tt':
            y3_AS.append(Obs_key['y3_AS'])
            y3_LO.append(Obs_key['y3_LO'])    
            wghts2_AS.append(Obs_key['wghts2_AS'])
            wghts2_LO.append(Obs_key['wghts2_LO'])        
            Q_AS.append(Obs_key['Q_AS'])
            Q_LO.append(Obs_key['Q_LO'])
            ydelta_AS.append(Obs_key['ydelta_AS'])
            ydelta_LO.append(Obs_key['ydelta_LO'])    

    '----=========================================================================----'
    '----=========================================================================----'
    '----==========================    HISTOGRAMMING    ==========================----'
    '----=========================================================================----'
    '----=========================================================================----'

    if FRAME == 'LAB':
        generate_normalized_histo(y3_AS, wghts_AS,'ACSy3', r'$y_{t}$', '$d$'r'$\sigma^{lab}_{A}$''/$d$' r'$y_{t}$ $[pb]$', (0, 2))
        generate_asymm_plot(y3_AS, y3_LO, wghts_AS, wghts_LO, 'FBy3', r'$y_{t}$', r'$A^{lab}_{FB}(y_{3})$', (0, 2))

    if FRAME == 'tt':
        
        generate_normalized_histo(Q_AS, wghts2_AS, 'ACSQ', r'$M_{t\bar{t}}$ $[GeV$ $]$', '$d$'r'$\sigma^{t\bar{t}}_{A}$''/$d$' r'$M_{t\bar{t}}$ $[fb/GeV$ $]$', (330, 800))
        generate_asymm_plot(Q_AS, Q_LO, wghts_AS, wghts_LO, 'FBQ', r'$M_{t\bar{t}}$ $[GeV$ $]$', r'$A^{lab}_{FB}$ $(M_{t\bar{t}})$', (330, 1200))
        
        generate_normalized_histo(ydelta_AS, wghts_AS, 'ACSydelta', r'$\Delta y$', '$d$'r'$\sigma^{t\bar{t}}_{A}$''/$d$' r'$\Delta y$ [pb]', (0,2))
        generate_asymm_plot(ydelta_AS, ydelta_LO, wghts_AS, wghts_LO, 'FBydelta', r'$\Delta y$', r'$A^{lab}_{FB}$ $(\Delta y)$', (0, 2))

