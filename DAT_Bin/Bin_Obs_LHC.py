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

#comment out the following if not using matplotlib and numpy
import matplotlib
import mpmath as mp
import numpy as np
import pylab as pl
import scipy
from scipy import interpolate, signal
import matplotlib.font_manager as fm
from matplotlib.font_manager import FontProperties
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
from collections import defaultdict
import resource
import gc
from scipy.interpolate import spline
from matplotlib.ticker import AutoMinorLocator

PDF1 = 'cteq66'
PDF2 = 'MSTW2008CPdeutnlo68cl'
PDF3 = 'NNPDF21_lostar_as_0130_100'
PDF_set = [PDF1, PDF2, PDF3]

 #Define variables in Kuhn paper
M_top = 173.1 # Mass of top quark GeV
alpha_QCD = 0.109 # alpha QCD
E_cut = 0.01 * M_top
deltath = 2

def mem():
    print('Memory usage         : % 8.3f GB' % (
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/10**6))
    
'----=========================================================================----'
'----==========================      LOAD FILES     ==========================----'
'----=========================================================================----'

print time.strftime("%H:%M")

cwd = os.getcwd()
base = os.path.splitext(os.path.basename(__file__))[0]

#get pickle names
pickle_names = []

for pickle_name in glob.glob('*.pickle'):
    dict_name = os.path.splitext(pickle_name)[0]
    pickle_names.append(dict_name)

def Tree(d, indent=0):
    '''Fucnction draws out dictionray tree'''
    for key, value in d.items():
      print('\t' * indent + str(key))
      if isinstance(value, dict):
         Tree(value, indent+1)
      else:
         print('\t' * (indent+1) + str(value)+'\n')

groupDict = {}

for pickle_key in pickle_names:
    #split key
    ss = pickle_key.split('_')
    mu_f = ss[0]+'_'+ss[1]
    FRAME = ss[2]+'_'+ss[3]
    CUT = ss[4]
    region = ss[5]
    
    if FRAME not in groupDict:
        groupDict[FRAME] = {}
        
    if mu_f not in groupDict[FRAME]:
        groupDict[FRAME][mu_f] = {}

    if CUT not in groupDict[FRAME][mu_f]:
        groupDict[FRAME][mu_f][CUT] = {}
        groupDict[FRAME][mu_f][CUT] = [pickle_key]
    else:
        groupDict[FRAME][mu_f][CUT].append(pickle_key)  

##Tree(groupDict)

'----=========================================================================----'
'----=========================================================================----'
'----=========================      OBSERVABLES     ==========================----'
'----=========================================================================----'
'----=========================================================================----'

FRAME_lA = ['ABS_y','PSE_RAP']
FRAME_lB = ['OUT_CUT', 'IN_CUT', 'OUT_CUTC']
FRAME_lC = ['Y_CUT',  'Y_CUTP', 'Y_CUTQ']

E_list = []
CUT_list = []

PLOT_d = {}

print 'Calculating asymmetries'
with open('C_ASYM.txt', 'w') as myfile:
    myfile.write('\n CHARGE ASYMMETRY ')
    
                
for FRAME_key in groupDict.keys():
    print '\n==============================================================================='
    print '===========================\t', FRAME_key,'\t==========================='
    print '===============================================================================\n'
    print time.strftime("%H:%M"), '\n'
    #create mu_key dict
    if FRAME_key not in PLOT_d:
        PLOT_d[FRAME_key] = {}
        
    for mu_key in groupDict[FRAME_key].keys():
        mem()
        print '\n=========================='
        print '===      ', mu_key,'     ==='
        print '==========================\n'
        #create FRAME_key dict
        if mu_key not in PLOT_d[FRAME_key]:
            PLOT_d[FRAME_key][mu_key] = {}

        Ccount = 0
        for CUT_key in groupDict[FRAME_key][mu_key].keys():
            Ccount +=1
            tot_Ccount = len(groupDict[FRAME_key][mu_key])
            
            #create LO CS
            if 'CUT' not in PLOT_d[FRAME_key][mu_key]:
                PLOT_d[FRAME_key][mu_key]['CUT'] = []
            #create LO CS
            if 'LO_CS' not in PLOT_d[FRAME_key][mu_key]:
                PLOT_d[FRAME_key][mu_key]['LO_CS'] = []
            #create AS CS            
            if 'AS_CS' not in PLOT_d[FRAME_key][mu_key]:
                PLOT_d[FRAME_key][mu_key]['AS_CS'] = []
            #create Charge asymmetry                            
            if 'C_Asymm' not in PLOT_d[FRAME_key][mu_key]:
                PLOT_d[FRAME_key][mu_key]['C_Asymm'] = []

            LO_CS = 0
            AS_CS = 0
            
            Pcount = 0
            for pickle_key in groupDict[FRAME_key][mu_key][CUT_key]:
                Pcount += 1
                tot_Pcount = len(groupDict[FRAME_key][mu_key][CUT_key])
                print ('loading Pickle file no. %s / %s for CUT: %s: %s / %s ...\n' % (Pcount, tot_Pcount, CUT_key, Ccount, tot_Ccount))

                '===================      OBSERVABLES     ======================'
                #load pickle file
                pickle_name = pickle_key+'.pickle'
                pickle_in = open(pickle_name, 'rb')
                data = pickle.load(pickle_in)
                print ('%s loaded!!\n' % pickle_key)

                if CUT_key not in CUT_list:
                    CUT_list.append(CUT_key)
                
                #define variables
                if FRAME_key in FRAME_lA:
                    CUT = CUT_key
                else:
                    CUT = float(CUT_key)
                    
                FRAME = FRAME_key
                int_region = pickle_key.split('_')[-1]
                mu_f = mu_key

                # Define ECM value 
                if FRAME in FRAME_lA:
                    #condition value for E_name
                    ECM = float(CUT_key[0]) * 1000
                    E_list.append(ECM)
                    s = ECM**2

                # For CUT need to define ECM 
                if FRAME in FRAME_lB:
                    ECM = 7000 
                    s = ECM**2

                # For CUT need to define ECM 
                if FRAME in FRAME_lC:
                    ECM = 7000
                    s = ECM**2 
                    
                tau0 = 4*M_top**2/s
                
                '----=========================================================================----'
                '----========================     qq CALCULATION    ==========================----'
                '----=========================================================================----'

                if int_region == 'qqbar':
                    qq_bar = data
                    print '----=========================================================================---- \n'
                    print ('Generating qq_bar observables, CUT = %s, mu = %s \n' % (CUT, mu_f))
                    #Random numbers from VEGAS
                    
                    wght = qq_bar.get('wghts')
                    LO_CS += sum(wght)

                '----=========================================================================----'
                '----========================     gg CALCULATION    ==========================----'
                '----=========================================================================----'  

                if int_region == 'gg':
                    gg = data
                    print '----=========================================================================---- \n'
                    print ('Generating gg observables, CUT = %s, mu = %s \n' % (CUT, mu_f))
                    #Random numbers from VEGAS
                    
                    wght = gg.get('wghts')
                    LO_CS += sum(wght)

                '----=========================================================================----'
                '----========================       SOFT REGION     ==========================----'
                '----=========================================================================----'  

                if int_region == 'soft':
                    soft = data
                    print '----=========================================================================---- \n'
                    print ('Generating soft observables, CUT = %s, mu = %s \n' % (CUT, mu_f))
                    #Random numbers from VEGAS
                    
                    wght = soft.get('wghts')
                    AS_CS += sum(wght)                

                '----=========================================================================----'
                '----========================       HARD REGION     ==========================----'
                '----=========================================================================----'  

                if int_region == 'hard':
                    hard = data
                    print '----=========================================================================---- \n'
                    print ('Generating hard observables, CUT = %s, mu = %s \n' % (CUT, mu_f))
                    #Random numbers from VEGAS
                    
                    wght = hard.get('wghts')
                    AS_CS += sum(wght)

                '----=========================================================================----'
                '----========================        QG REGION      ==========================----'
                '----=========================================================================----'  

                if int_region == 'qg':
                    qg = data
                    print '----=========================================================================---- \n'
                    print ('Generating qg observables, CUT = %s, mu = %s \n' % (CUT, mu_f))
                    #Random numbers from VEGAS
                    
                    wght = qg.get('wghts')
                    AS_CS += sum(wght)

                pickle_in.close()
                mem()
                NEVAL = len(wght)
                
            
            '----===========================    CHARGE ASYMMETRY    ==========================----'
            
            print 'Calculating asymmetries'
            with open('C_ASYM.txt', 'a') as myfile:
                myfile.write('\n CHARGE ASYMMETRY (%s frame, CUT = %s)' % (FRAME, CUT))

            if LO_CS ==0:
                C_Asymm = 0
            else:
                C_Asymm = AS_CS/LO_CS * 100
            print '\n=============================='
            print ' || A_CS: %s' % AS_CS
            print ' || LO_CS: %s' % LO_CS
            print ' || C Asymmetry: %s %%' % C_Asymm
            print ' || NEVAL: %s' % NEVAL
            print ' || mu_f: %s' % mu_f
            print ' || CUT: %s' % CUT            
            print '==============================\n'

            with open('C_ASYM.txt', 'a') as myfile:
                myfile.write('\n\n============================\n')
                myfile.write('\n CUT: %s || mu_f: %s\n\n' % (CUT, mu_key))
                myfile.write('|| A_CS: %s\n|| LO_CS: %s\n|| C Asymmetry: %s %%\n|| NEVAL: %s\n' % (AS_CS, LO_CS, C_Asymm, NEVAL))
                myfile.write('============================\n')
                
            add = {'AS_CS': AS_CS, 'LO_CS': LO_CS, 'C_Asymm': C_Asymm, 'CUT': CUT}

            for key, val in add.items():
                if key in PLOT_d[FRAME_key][mu_key]:
                    PLOT_d[FRAME_key][mu_key][key].append(val)
                else:
                   'Error!!! ERROR!!! ERRRRRROOOOOOOOORRRR!!!!'
                            
##Tree(PLOT_d)

'----=========================================================================----'
'----==========================      SAVE GRAPH     ==========================----'
'----=========================================================================----'
 #make directory path

Save_path = cwd+'/GRAPHS/'

 #makes directories
if not os.path.exists(os.path.dirname(Save_path)):
    try:
        os.makedirs(os.path.dirname(Save_path))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise
        
'----=========================================================================----'
'----========================      CREATE TABLE      =========================----'
'----=========================================================================----'
GRAPH = {}

for FRAME_key in PLOT_d:
    CUT = PLOT_d[FRAME_key]['M_topp'].get('CUT')
    AS_CS = PLOT_d[FRAME_key]['M_topp'].get('AS_CS')
    LO_CS = PLOT_d[FRAME_key]['M_topp'].get('LO_CS')
    C_Asymm = PLOT_d[FRAME_key]['M_topp'].get('C_Asymm')

    #sort in order
    if FRAME_key in FRAME_lA:
        print CUT, '0'
        CUT = [float(x[:-3]) for x in CUT]
        print CUT, '1st'
        CUT, AS_CS, LO_CS, C_Asymm = zip(*sorted(zip(CUT, AS_CS, LO_CS, C_Asymm)))
        print CUT, '2nd'
        CUT = [str(x)[:-2]+'TeV' for x in CUT]
        print CUT


    else: 
        CUT, AS_CS, LO_CS, C_Asymm = zip(*sorted(zip(CUT, AS_CS, LO_CS, C_Asymm)))
    
    if C_Asymm == None:
        print 'M_topp values not found ...'
        break

    #create array for M_half
    if PLOT_d[FRAME_key].get('M_half') != None:
        ref = PLOT_d[FRAME_key]['M_half']

        CUT_up = ref['CUT']
        AS_CS_up = ref['AS_CS']
        LO_CS_up = ref['LO_CS']
        C_Asymm_up = ref['C_Asymm']

        CUT_up, AS_CS_up, LO_CS_up, C_Asymm_up = zip(*sorted(zip(CUT_up, AS_CS_up, LO_CS_up, C_Asymm_up)))

        v0 = zip(CUT, AS_CS, LO_CS, C_Asymm)
        vup = zip(CUT_up, AS_CS_up, LO_CS_up, C_Asymm_up)       

        #calculate upper errors
        up_E = []
        for i in range(len(v0)):
            if v0[i][0] in CUT_up:
                index = CUT_up.index(v0[i][0])
                up_E.append([x - y for x, y in zip(vup[index][1:], v0[i][1:])])
            else:
                up_E.append([v0[i][0], None, None, None])
                
        AS_CS_up_E, LO_CS_up_E, C_Asymm_up_E = zip(*up_E)

    #create array for M_doub
    if PLOT_d[FRAME_key].get('M_doub') != None:
        ref = PLOT_d[FRAME_key]['M_doub']

        CUT_low = ref['CUT']
        AS_CS_low = ref['AS_CS']
        LO_CS_low = ref['LO_CS']
        C_Asymm_low = ref['C_Asymm']

        v0 = zip(CUT, AS_CS, LO_CS, C_Asymm)
        vup = zip(CUT_low, AS_CS_low, LO_CS_low, C_Asymm_low)       

        #calculate lower errors
        low_E = []        
        for i in range(len(v0)):
            if v0[i][0] in CUT_low:
                index = CUT_low.index(v0[i][0])
                low_E.append([x - y for x, y in zip(vup[index][1:], v0[i][1:])])
            else:
                low_E.append([v0[i][0], None, None, None])
                
        AS_CS_low_E, LO_CS_low_E, C_Asymm_low_E = zip(*low_E)        

        '----===========================    DRAW TABLE    ==========================----'

        if FRAME == 'ABS_y':
##            rows = [r'$LHC %s TeV$' % x for x in CUT]
            columns = ('', r'$\sigma_{A}^{lab}$', r'$\sigma_{S}^{lab}$', r'$A_{C}^{y} (\%)$')

        if FRAME == 'PSE_RAP':
##            rows = [r'$LHC %s TeV$' % x for x in CUT]
            columns = ('', r'$\sigma_{A}^{lab}$', r'$\sigma_{S}^{lab}$', r'$A_{C}^{\nu} (\%)$')
            
        elif FRAME in FRAME_lC:
##            rows = [r'$%s$' % x for x in CUT]
            columns = (r'$Y_{cut}$', r'$\sigma_{A}^{lab}$' '$[pb]$', r'$\sigma_{S}^{lab}$' '$[pb]$', r'$A_{C}^{t\bar{t}} (\%)$')
            
        elif FRAME == 'OUT_CUT' or FRAME == 'OUT_CUTC':
##            rows = [r'$%s$' % x for x in CUT]
            columns = (r'$y_{C}^{out}$', r'$\sigma_{A}^{lab}$' '$[pb]$', r'$\sigma_{S}^{lab}$' '$[pb]$', r'$A_{C}^{out} (\%)$')

        elif FRAME == 'IN_CUT':
##            rows = [r'$%s$' % x for x in CUT]
            columns = (r'$y_{C}^{in}$', r'$\sigma_{A}^{lab}$' '$[pb]$', r'$\sigma_{S}^{lab}$' '$[pb]$', r'$A_{C}^{in} (\%)$')


        fig, ax = plt.subplots()

        AS_zip = zip(AS_CS, AS_CS_up_E, AS_CS_low_E)
        LO_zip = zip(LO_CS, LO_CS_up_E, LO_CS_low_E)
        C_Asymm_zip = zip(C_Asymm, C_Asymm_up_E, C_Asymm_low_E)

        AS_zip = [r'$%.3f^{+%.3f}_{-%.3f}$' % x for x in AS_zip]
        LO_zip = [r'$%.3f^{+%.3f}_{-%.3f}$' % x for x in LO_zip]
        C_Asymm_zip = [r'$%.2f^{+%.2f}_{-%.2f}$' % x for x in C_Asymm_zip]        
        
        zipall = zip(AS_zip, LO_zip, C_Asymm_zip)

        table_values = []
        
        for i in range(len(zipall)):
            if FRAME_key in FRAME_lA:
                table_values.append(np.append(r'$%s$' % CUT[i][:-3]+ r' $TeV$', zipall[i]))
            else:
                table_values.append(np.append(r'$%.2f$' % CUT[i], zipall[i]))
        
        lightgrn = (0.5, 0.8, 0.5) #color of columns
        colors = [[(0.99, 0.99, 0.99) for c in range(len(columns))] for r in range(len(table_values))] #color of cells
                  
        ax.axis('tight')
        ax.axis('off')

        the_table =ax.table(cellText=table_values, cellColours=colors,
          cellLoc='center', colWidths=[0.2 for x in columns],
          rowLabels=None, rowColours=None, rowLoc='center',
          colLabels=columns, colColours=[lightgrn]*len(columns), colLoc='center',
          loc='center', bbox=None)
        
        the_table.scale(1, 2)
        fig.tight_layout()
        
        for key, cell in the_table.get_celld().items():
            row, col = key
            if row == 0:
                cell.set_text_props(fontproperties=FontProperties(weight='extra bold',style='italic', size=20))
            else:
                cell.set_text_props(fontproperties=FontProperties(weight='extra bold',style='italic', size=20))

        plt.savefig(Save_path + 'TABLE_' + FRAME_key + '.pdf')
        plt.show()
        
'----=========================================================================----'
'----=========================================================================----'
'----==========================       GRAPHING      ==========================----'
'----=========================================================================----'
'----=========================================================================----'


'----=========================================================================----'
'----======================        OUTCUT A CS      ==========================----'
'----=========================================================================----'
#Get OUT_CUT plotting variables
if PLOT_d.get('OUT_CUT') != None:
    
    CUT_axisO = np.asarray(PLOT_d['OUT_CUT']['M_topp'].get('CUT'))
    print CUT_axisO
    AS_CSO = np.asarray(PLOT_d['OUT_CUT']['M_topp'].get('AS_CS'))
    C_AsymmO = np.asarray(PLOT_d['OUT_CUT']['M_topp'].get('C_Asymm'))

    CUT_axisO, AS_CSO, C_AsymmO = zip(*sorted(zip(CUT_axisO, AS_CSO, C_AsymmO)))
    
    if CUT_axisO == None or AS_CSO == None:
        print 'error! Error!! ERROR!!!'

    
    if PLOT_d['OUT_CUT'].get('M_half') != None:
        CUT_axis_upO = np.asarray(PLOT_d['OUT_CUT']['M_half'].get('CUT'))
        AS_CS_upO = np.asarray(PLOT_d['OUT_CUT']['M_half'].get('AS_CS'))
        C_Asymm_upO = np.asarray(PLOT_d['OUT_CUT']['M_half'].get('C_Asymm'))

        CUT_axis_upO, AS_CS_upO, C_Asymm_upO = zip(*sorted(zip(CUT_axis_upO, AS_CS_upO, C_Asymm_upO)))


    if PLOT_d['OUT_CUT'].get('M_doub') != None:
        CUT_axis_lowO = np.asarray(PLOT_d['OUT_CUT']['M_doub'].get('CUT'))
        AS_CS_lowO = np.asarray(PLOT_d['OUT_CUT']['M_doub'].get('AS_CS'))
        C_Asymm_lowO = np.asarray(PLOT_d['OUT_CUT']['M_doub'].get('C_Asymm'))

        CUT_axis_lowO, AS_CS_lowO, C_Asymm_lowO = zip(*sorted(zip(CUT_axis_lowO, AS_CS_lowO, C_Asymm_lowO)))


    # Plot AS_CS figure
    fig, ax = plt.subplots()
    CUT_axisO = np.asarray(CUT_axisO)
    #smooth out curve
    xnew = np.linspace(np.min(CUT_axisO), np.max(CUT_axisO), 300)
    AS_CSO = spline(CUT_axisO, AS_CSO, xnew)
    
    if PLOT_d['OUT_CUT'].get('M_half') and PLOT_d['OUT_CUT'].get('M_doub') != None:
        
        AS_CS_upO = np.asarray(AS_CS_upO)
        AS_CS_lowO = np.asarray(AS_CS_lowO)
        #smooth out curve
        AS_CS_upO = spline(CUT_axisO, AS_CS_upO, xnew)
        AS_CS_lowO = spline(CUT_axisO, AS_CS_lowO, xnew)

        plt.plot(xnew, AS_CS_upO, label = r'$M_{top}/$ $2$', color='black', linewidth=2.5, alpha=1, linestyle= '-')
        plt.plot(xnew, AS_CSO, label = r'$M_{top}$', color='black', linewidth=1.6, alpha=1, linestyle= '--')
        plt.plot(xnew, AS_CS_lowO, label = r'$M_{top}$'r'$\times$ $2$', color='black', linewidth=0.9, alpha=1, linestyle= '-')

        ax.fill_between(xnew, AS_CS_upO, AS_CS_lowO, where=AS_CS_upO >= AS_CS_lowO, facecolor='grey', interpolate=True)

    else:
        plt.plot(xnew, AS_CSO, label = r'$M_{top}$', color='black', linewidth=1.6, alpha=1, linestyle= '-')    

    #set minor axis
    minorLocator = AutoMinorLocator(5)
    ax.xaxis.set_minor_locator(minorLocator)
    minorLocator = AutoMinorLocator(5)
    ax.yaxis.set_minor_locator(minorLocator)

    # update the view limits
    ax.set_xlim(1, 2.5)
##    ax.set_ylim(0, 2.5)
 
    
    plt.legend(frameon=False, fontsize = 16)
    plt.xlabel(r'$y_{cut}$', fontsize=20)
    plt.ylabel(r'$\sigma_{A}^{lab}$' r'$(y_{cut})$', fontsize=22)
    plt.tight_layout()

    plt.savefig(Save_path + 'A_CS_' + FRAME_key + '.pdf')
##    plt.show()


'----=========================================================================----'
'----======================        C ASYMM PLOT      =========================----'
'----=========================================================================----'

if PLOT_d.get('IN_CUT') or PLOT_d.get('OUT_CUT') or PLOT_d.get('OUT_CUTC')!= None:
    
    fig, ax = plt.subplots()    

    #Get IN_CUT plotting variables
    if PLOT_d.get('IN_CUT') != None:
        
        CUT_axisI = np.asarray(PLOT_d['IN_CUT']['M_topp'].get('CUT'))
        AS_CSI = np.asarray(PLOT_d['IN_CUT']['M_topp'].get('AS_CS'))
        C_AsymmI = np.asarray(PLOT_d['IN_CUT']['M_topp'].get('C_Asymm'))

        CUT_axisI, AS_CSI, C_AsymmI = zip(*sorted(zip(CUT_axisI, AS_CSI, C_AsymmI)))


        if CUT_axisI == None or AS_CSI == None:
            print 'error! Error!! ERROR!!!'

        
        if PLOT_d['IN_CUT'].get('M_half') != None:
            CUT_axis_upI = np.asarray(PLOT_d['IN_CUT']['M_half'].get('CUT'))
            AS_CS_upI = np.asarray(PLOT_d['IN_CUT']['M_half'].get('AS_CS'))
            C_Asymm_upI = np.asarray(PLOT_d['IN_CUT']['M_half'].get('C_Asymm'))

            CUT_axis_upI, AS_CS_upI, C_Asymm_upI = zip(*sorted(zip(CUT_axis_upI, AS_CS_upI, C_Asymm_upI)))


        if PLOT_d['IN_CUT'].get('M_doub') != None:
            CUT_axis_lowI = np.asarray(PLOT_d['IN_CUT']['M_doub'].get('CUT'))
            AS_CS_lowI = np.asarray(PLOT_d['IN_CUT']['M_doub'].get('AS_CS'))
            C_Asymm_lowI = np.asarray(PLOT_d['IN_CUT']['M_doub'].get('C_Asymm'))

            CUT_axis_lowI, AS_CS_lowI, C_Asymm_lowI = zip(*sorted(zip(CUT_axis_lowI, AS_CS_lowI, C_Asymm_lowI)))
            

    #Get IN_CUT plotting variables
    if PLOT_d.get('OUT_CUTC') != None:
        
        CUT_axisC = np.asarray(PLOT_d['OUT_CUTC']['M_topp'].get('CUT'))
        AS_CSC = np.asarray(PLOT_d['OUT_CUTC']['M_topp'].get('AS_CS'))
        C_AsymmC = np.asarray(PLOT_d['OUT_CUTC']['M_topp'].get('C_Asymm'))

        CUT_axisC, AS_CSC, C_AsymmC = zip(*sorted(zip(CUT_axisC, AS_CSC, C_AsymmC)))


        if CUT_axisC == None or AS_CSC == None:
            print 'error! Error!! ERROR!!!'

        
        if PLOT_d['OUT_CUTC'].get('M_half') != None:
            CUT_axis_upC = np.asarray(PLOT_d['OUT_CUTC']['M_half'].get('CUT'))
            AS_CS_upC = np.asarray(PLOT_d['OUT_CUTC']['M_half'].get('AS_CS'))
            C_Asymm_upC = np.asarray(PLOT_d['OUT_CUTC']['M_half'].get('C_Asymm'))

            CUT_axis_upC, AS_CS_upC, C_Asymm_upC = zip(*sorted(zip(CUT_axis_upC, AS_CS_upC, C_Asymm_upC)))
            

        if PLOT_d['OUT_CUTC'].get('M_doub') != None:
            CUT_axis_lowC = np.asarray(PLOT_d['OUT_CUTC']['M_doub'].get('CUT'))
            AS_CS_lowC = np.asarray(PLOT_d['OUT_CUTC']['M_doub'].get('AS_CS'))
            C_Asymm_lowC = np.asarray(PLOT_d['OUT_CUTC']['M_doub'].get('C_Asymm'))

            CUT_axis_lowC, AS_CS_lowC, C_Asymm_lowC = zip(*sorted(zip(CUT_axis_lowC, AS_CS_lowC, C_Asymm_lowC)))


    fcolour_list = []
    labelit_list = []

    #Plot IN, OUT, MCUT
    if PLOT_d.get('OUT_CUT') != None:
        
        fcolour_list.append('red')
        labelit_list.append(r'$A_{C}^{out}$')
        
        CUT_axisO = np.asarray(CUT_axisO)
        C_AsymmO = np.asarray(C_AsymmO)
        C_Asymm_upO = np.asarray(C_Asymm_upO)
        C_Asymm_lowO = np.asarray(C_Asymm_lowO)

        #smooth out curve
        xnew = np.linspace(np.min(CUT_axisO), np.max(CUT_axisO), 300)
        C_AsymmO = spline(CUT_axisO, C_AsymmO, xnew)
        C_Asymm_upO = spline(CUT_axisO, C_Asymm_upO, xnew)
        C_Asymm_lowO = spline(CUT_axisO, C_Asymm_lowO, xnew)
        
        plt.plot(xnew, C_AsymmO, color='black', linewidth=1.6, alpha=1, linestyle= '--')
        plt.plot(xnew, C_Asymm_upO, color='black', linewidth=1.6, alpha=1, linestyle= '-')
        plt.plot(xnew, C_Asymm_lowO, color='black', linewidth=1.6, alpha=1, linestyle= '-')

        ax.fill_between(xnew, C_Asymm_upO, C_Asymm_lowO, where=C_Asymm_upO >= C_Asymm_lowO, facecolor='red', interpolate=True)    


    if PLOT_d.get('IN_CUT') != None:
        
        fcolour_list.append('blue')
        labelit_list.append(r'$A_{C}^{in}$')
        
        CUT_axisI = np.asarray(CUT_axisI)
        C_AsymmI = np.asarray(C_AsymmI)        
        C_Asymm_upI = np.asarray(C_Asymm_upI)
        C_Asymm_lowI = np.asarray(C_Asymm_lowI)

        #smooth out curve
        xnew = np.linspace(np.min(CUT_axisI), np.max(CUT_axisI), 300)
        C_AsymmI = spline(CUT_axisI, C_AsymmI, xnew)
        C_Asymm_upI = spline(CUT_axisI, C_Asymm_upI, xnew)
        C_Asymm_lowI = spline(CUT_axisI, C_Asymm_lowI, xnew)
        
        plt.plot(xnew, C_AsymmI, color='black', linewidth=1.6, alpha=1, linestyle= '--')
        plt.plot(xnew, C_Asymm_upI, color='black', linewidth=1.6, alpha=1, linestyle= '-')
        plt.plot(xnew, C_Asymm_lowI, color='black', linewidth=1.6, alpha=1, linestyle= '-')

        ax.fill_between(xnew, C_Asymm_upI, C_Asymm_lowI, where=C_Asymm_upI >= C_Asymm_lowI, facecolor='blue', interpolate=True)    


    if PLOT_d.get('OUT_CUTC') != None:

        fcolour_list.append('green')
        labelit_list.append(r'$A_{C}^{out}$' r'$(|y_{t}|,$' r'$|y_{\bar{t}}|$' r'$<\ 2.5$')
        
        CUT_axisC = np.asarray(CUT_axisC)
        C_AsymmC = np.asarray(C_AsymmC)    
        C_Asymm_upC = np.asarray(C_Asymm_upC)
        C_Asymm_lowC = np.asarray(C_Asymm_lowC)

        #smooth out curve
        xnew = np.linspace(np.min(CUT_axisC), np.max(CUT_axisC), 300)
        C_AsymmC = spline(CUT_axisC, C_AsymmC, xnew)
        C_Asymm_upC = spline(CUT_axisC, C_Asymm_upC, xnew)
        C_Asymm_lowC = spline(CUT_axisC, C_Asymm_lowC, xnew)    
        
        plt.plot(xnew, C_AsymmC, color='black', linewidth=1.6, alpha=1, linestyle= '--')
        plt.plot(xnew, C_Asymm_upC, color='black', linewidth=1.6, alpha=1, linestyle= '-')
        plt.plot(xnew, C_Asymm_lowC, color='black', linewidth=1.6, alpha=1, linestyle= '-')

        ax.fill_between(xnew, C_Asymm_upC, C_Asymm_lowC, where=C_Asymm_upC >= C_Asymm_lowC, facecolor='green', interpolate=True)    


     #draw legend
    legend_patch = [patches.Patch(color=fcolour_list[i], label = "{:s}".format(labelit_list[i])) for i in range(len(labelit_list))]
    plt.legend(handles=legend_patch, bbox_to_anchor=(0.5, 0.97), ncol=1, frameon=False, fontsize = 16)

     #set font dict
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 16,
            }


    #set minor axis
    minorLocator = AutoMinorLocator(5)
    ax.xaxis.set_minor_locator(minorLocator)
    minorLocator = AutoMinorLocator(5)
    ax.yaxis.set_minor_locator(minorLocator)

    plt.xlabel(r'$y_{cut}$', fontsize=20)
    plt.ylabel(r'$A_{C}^{in/out}$' r'$(\%)(y_{cut})$', fontsize=22)

    #print LHC power
    ##xp = array_range[0] + (array_range[1] - array_range[0]) / 3
    ##yp = np.max(top_lim) - np.max(top_lim)/20
    ##plt.text(xp, yp, r'$\sqrt{S}=1.96\ TeV$' , fontdict=font) #r'$math.sqrt(S)=1.96\ TeV$'

    plt.savefig(Save_path + 'C_Asymm_' + FRAME_key + '.pdf')

    plt.tight_layout()
    plt.show()

        
    
