#!/usr/bin/env python
#
# AUTHOR:
# Rostam Golesorkhtabar 
# r.golesorkhtabar@gmail.com
# 
# DATE:
# Tue Jan 01 00:00:00 2013
#
# EXPLANATION:
# The parameters in equation of states:
# --- Murnaghan ---
#     E(V) = E0 + B0*V/B' * [(V0/V)^B' / (B'-1) + 1] - B0 V0 / (B'-1)
#     P(V) = B0/B'*((V0/V)**B' -1)
#
# --- Birch-Murnaghan ---
#     E(V) = E0 + 9B0V0/16{[(V0/V)^(2/3)-1]^3 B' + [(V0/V)^(2/3)-1]^2 [6-4(V0/V)^(2/3)]}
#     P(V) = 3/2*B0*((V0/V)**(7/3) - (V0/V)**(5/3))*(1 + 3/4*(B'-4)*[(V0/V)**(2/3) - 1])
#
# E0: Energy in equilibrium
# V0: Volume in equilibrium
# B0: Bulk modulus in equilibrium
# B': Pressure derivative of B
#
# See also: http://en.wikipedia.org/wiki/Birch-Murnaghan_equation_of_state
#__________________________________________________________________________________________________

from   pylab import *
import os
import sys
import glob
import copy
import math
import os.path
import numpy as np
from   lxml  import etree as ET
import matplotlib.pyplot as plt
import pylab             as pyl
from   scipy.optimize import fmin_powell

#%!%!%--- CONSTANTS ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
_e    =  1.602176565e-19          # elementary charge
Bohr  =  5.291772086e-11          # Bohr to meter
Ryd2eV= 13.605698066              # Ryd to eV
Ha2Ryd=  2.                       # Ha to Ryd
ToGPa = (_e*Ryd2eV)/(1e9*Bohr**3) # Ryd/[Bohr]^3 to GPa
#__________________________________________________________________________________________________

#%!%!%--- SUBROUTINS AND FUNCTIONS ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
def E_eos(p0, V):
    if (eos=='M'):
        """ Murnaghan Energy"""
        E0, V0, B0, Bp = p0
        E = E0 + (B0*V/Bp*(1/(Bp-1)*(V0/V)**Bp +1)-B0*V0/(Bp-1))
    else:
        """ Birch-Murnaghan Energy"""
        E0, V0, B0, Bp = p0
        E = E0 + (9.*B0*V0/16)*(((((V0/V)**(2./3))-1.)**3.)*Bp \
               + ((((V0/V)**(2/3.))-1.)**2.)*(6.-4.*((V0/V)**(2./3.))))
    return E
#--------------------------------------------------------------------------------------------------

def P_eos(p0, V):
    if (eos=='M'):
        """ Murnaghan Pressure"""
        E0, V0, B0, Bp = p0
        P = B0/Bp*((V0/V)**Bp - 1.)
    else:
        """ Birch-Murnaghan Pressure"""
        E0, V0, B0, Bp = p0
        P = 3./2*B0*((V0/V)**(7./3) - (V0/V)**(5./3))*(1. + 3./4*(Bp-4.)*((V0/V)**(2./3) - 1.))
    return P
#--------------------------------------------------------------------------------------------------

def snr(p0, v, e):
    """ Squared norm of residue vector calculation """
    return np.sum((e - E_eos(p0, v))**2.)
#--------------------------------------------------------------------------------------------------

def sortlist(lst1, lst2):
    temp = copy.copy(lst1)

    lst3 = []
    lst4 = []

    temp.sort()

    for i in range(len(lst1)):
        lst3.append(lst1[lst1.index(temp[i])])
        lst4.append(lst2[lst1.index(temp[i])])

    return lst3, lst4
#__________________________________________________________________________________________________

#%!%!%--- READING SECTION ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
print'                                                                    \
\n     +-----------------------------------------------------------------+\
\n     |*****************************************************************|\
\n     |*                                                               *|\
\n     |*                WELCOME TO THE "eos.py" PROGRAM                *|\
\n     |*       Rostam Golesorkhtabar, r.golesorkhtabar@gmail.com       *|\
\n     |*                                                               *|\
\n     |*****************************************************************|\
\n     +-----------------------------------------------------------------+'

INF = raw_input('\n>>>> Enter "volume energy" file name: ')
if (os.path.exists(INF)==False):
    sys.exit('\n    ... Oops ERROR: There is NO '+ INF +' file !?!?!?    \n')

v_array, e_array = np.loadtxt(INF).T
v_list , e_list  = sortlist(list(v_array), list(e_array))
vi = array(v_list)
ei = array(e_list)
if (len(ei) < 3): sys.exit('\n    ... Oops ERROR: EOS fit needs at least 3 points.    \n')

print '\
\n     Please specify the "volume unit".\
\n     [Bohr^3] ---=> 1                 \
\n     [Ang.^3] ---=> 2                 '
num = input(">>>> Choose '1' or '2': ")
if (num != 1 and num != 2 ): sys.exit("\n    ...Oops ERROR: Choose '1' or '2' \n")
if (num == 1): v_unit = 'Bohr^3'; v2Bhr= 1.
if (num == 2): v_unit = 'Ang.^3'; v2Bhr= (Bohr*1e+10)**-3.

print '\
\n     Please specify the "energy unit".\
\n     Rydberg  -------=> 1             \
\n     Hartree  -------=> 2             \
\n     electron Volt --=> 3             '
num = input(">>>> Choose '1', '2', or '3': ")
if (num != 1 and num != 2 and num != 3 ):
    sys.exit("\n    ...Oops ERROR: Choose '1', '2', or '3'  \n")
if (num == 1): e_unit = 'Ryd'; e2Ryd= 1.
if (num == 2): e_unit = 'Ha' ; e2Ryd= Ha2Ryd
if (num == 3): e_unit = 'eV' ; e2Ryd= 1./Ryd2eV

vi = vi*v2Bhr
ei = ei*e2Ryd

#%!%!%!%!--- Reading the EOS type ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
eos = raw_input('\n>>>> Murnaghan or Birch-Murnaghan EOS: [M/B] ').upper()
if (eos != 'B' and eos != 'M'): sys.exit("\n    ... Oops ERROR: Choose 'B' or 'M' \n")
if (eos == 'B'): eos = 'BM'
#--------------------------------------------------------------------------------------------------

#%!%!%--- FIT CALCULATIONS ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
a2, a1, a0 = np.polyfit(vi, ei, 2)
V0 = -a1/(2.*a2)
E0 = a2*V0**2. + a1*V0 + a0
B0 = a2*V0
Bp = 2.0
p0 = [E0, V0, B0, Bp]

viei = sorted([zip(vi, ei)])
v, e = np.array(viei).T

print
p1, fopt, direc, n_iter, n_funcalls, warnflag = \
fmin_powell(snr, p0, args=(v, e), full_output=True)
E0, V0, B0, Bp = p1

E0_dim = E0/e2Ryd
V0_dim = V0/v2Bhr
B0_GPa = B0*ToGPa
fopt = fopt/(e2Ryd**2.)

print\
'\n Log(Final Residue in '+e_unit+'):'+str(round(log10(sqrt(fopt)),3)),'\n'\
'\n === Final Parameters ===='                 \
'\n E0 = '+str(round(E0_dim,8))+' ['+e_unit+']'\
'\n V0 = '+str(round(V0_dim,4))+' ['+v_unit+']'\
'\n B0 = '+str(round(B0_GPa,3))+' [GPa]'       \
"\n B' = "+str(round(Bp,3)),                   \
'\n ========================='
#--------------------------------------------------------------------------------------------------

#%!%!%--- WRITING THE OUTPUT FILE ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
outf = open(eos+'_EOS.out', 'w')

if (eos=='M'):
    print >>outf,' === Murnaghan ==========='
else:
    print >>outf,' === Birch-Murnaghan ====='

print >>outf,                                  \
  ' E0 = '+str(round(E0_dim,8))+' ['+e_unit+']'\
'\n V0 = '+str(round(V0_dim,4))+' ['+v_unit+']'\
'\n B0 = '+str(round(B0_GPa,3))+' [GPa]'       \
"\n B' = "+str(round(Bp,3))    +               \
'\n ========================='                 \
'\n Volume    E_dft-E_eos    P_eos [GPa]    -dE/dV [GPa]'

for i in range(len(ei)):
    Pi_eos = P_eos(p1, vi[i])*ToGPa
    ei_eos = E_eos(p1, vi[i])

    if (Pi_eos > 0):
        Peos = '+'+str(round(Pi_eos,3))
    else:
        Peos =     str(round(Pi_eos,3))

    if (i+1 < len(ei)): 
        Pi_dft = -(ei[i+1]-ei[i])*ToGPa/(vi[i+1]-vi[i])

    if (Pi_dft > 0):
        Pdft = '+'+str(round(Pi_dft,3))
    else:
        Pdft =     str(round(Pi_dft,3))

    if (i+1 == len(ei)): Pdft = ' NaN'

    print >>outf, str(round(vi[i]/v2Bhr,4))       \
                , '%12.8f'%((ei[i]-ei_eos)/e2Ryd) \
                , '    '                          \
                , Peos                            \
                , '    '                          \
                , Pdft
outf.close()
#--------------------------------------------------------------------------------------------------

#%!%!%--- PLOT SECTION ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
vmn  = min(np.min(vi), V0)
vmx  = max(np.max(vi), V0)
dv   = vmx - vmn
v_eos= np.linspace(vmn-(0.1*dv), vmx+(0.1*dv), 1000)
e_eos= E_eos(p1, v_eos)/e2Ryd
v_eos= v_eos/v2Bhr
ei   = ei/e2Ryd
vi   = vi/v2Bhr

fig = plt.figure()
ax  = fig.add_subplot(111)
fig.subplots_adjust(left=0.20)

plt.text(0.3,0.75,'E$_{min}$ = '  +str(round(E0_dim,8))+' ['+e_unit+']',transform = ax.transAxes)
plt.text(0.3,0.70,'V$_{min}$ = '  +str(round(V0_dim,4))+' ['+v_unit+']',transform = ax.transAxes)
plt.text(0.3,0.65,'B$_0$ = '      +str(round(B0_GPa,3))+' [GPa]'       ,transform = ax.transAxes)
plt.text(0.3,0.60,'B$^\prime$ = '+str(round(Bp,3))                    ,transform = ax.transAxes)

xx = [] ; xx = v_eos
yy = [] ; yy = e_eos
x0 = [] ; x0 = vi
y0 = [] ; y0 = ei

if (eos=='M'):
    eos_label='Murnaghan EOS'
else:
    eos_label='Birch-Murnaghan EOS'

ax.set_xlabel('Volume ['+v_unit+']', fontsize = 18)
ax.set_ylabel('Energy ['+e_unit+']', fontsize = 18)

ax.plot(xx, yy, 'k'               ,
                color     = 'red' ,
                linewidth = 2     ,
                label     = eos_label)

ax.plot(x0, y0, 'o'                     ,
                color          = 'green',
                markersize     = 8      ,
                markeredgecolor= 'black',
                markeredgewidth= 1      ,
                label          = 'DFT Calc.')
ax.legend(numpoints=1,loc=9)

for label in ax.xaxis.get_ticklabels(): label.set_fontsize(15)
for label in ax.yaxis.get_ticklabels(): label.set_fontsize(15)
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(6)
    line.set_markeredgewidth(2)

pyl.grid(True)   
ax.xaxis.set_major_locator(MaxNLocator(7))

dyy = (max(yy)-min(yy))/15
ax.set_ylim(min(yy)-dyy,max(yy)+dyy)
dxx = (max(xx)-min(xx))/18
ax.set_xlim(min(xx)-dxx,max(xx)+dxx)

plt.savefig(eos+'_EOS.png',orientation='portrait',format='png',dpi=300)
plt.savefig(eos+'_EOS.eps',orientation='portrait',format='eps')
plt.show()
#--------------------------------------------------------------------------------------------------
