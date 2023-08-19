#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 15:09:07 2021

@author: arventh
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

#Inputfile wtih time, energy, HBs, strand length and x, y, z co-ords for 2 the two traps, i.e. Zf[4:7], Zf[7:10]
folder = '/oxDNAdata/Rong/24mer_peeling' #raw_input directory
txt = '7-2_out_Observables_1.dat'
os.chdir(folder)
data = open(txt, 'r')
break_force, i, x, y, y_peak, err = 0,0,[],[],[], []

#Constants. Change these parameters to simulation conditons.
dt = 0.005 #simulation time units (Output has time in MD units where MD units = steps * dt)
F_units = 48.63 #units in pN 
l_units = 0.8518 #units in nm 
t_units = 3.03e-12 #seconds
k1,k2 = 0.2,0.2 #simulation units
m_avg = 4000 #exponential moving average smoothing of n-points
ext_rate = 0.1e-7 #rate of increase length (simulation units) per simulation step
F_dir = np.array([0,0,-1]) #direction of trap movement
t1_init = np.array([0,0,0]) #moving trap
t2_init = np.array([0,0,0]) #fixed trap
UnitF_dir = F_dir/np.linalg.norm(F_dir) #unit vector
init_dist = np.dot((t1_init - t2_init), UnitF_dir) #initial distance between two traps along force dir
keff = ((k1*k2)/(k1+k2)) #1 unit of force constant (1 unit force/1 unit length) - 57.09 pN/nm

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for line in data:
    Zf = np.array(line.split(), dtype=float)
    l_dna = Zf[3] #end to end strand length 
    trap1_current = t1_init + (Zf[0]/dt)*ext_rate*UnitF_dir #trap1 displacement
    trap_ext1 = trap1_current - Zf[7:10] #vector subtraction of moving trap1 and attached nucleotide
    trap_ext2 = Zf[4:7] - t2_init #vector subtraction of fixed trap2 and attached nucleotide
    force = (np.dot((trap_ext1+trap_ext2),UnitF_dir))*keff*F_units
    overall_ext = trap1_current - t2_init
    trap_ext = (np.dot(overall_ext,UnitF_dir)-init_dist)*l_units
    x.append(force) #this is extension of traps along force-axis not the DNA strand extension.
    y.append(Zf[2])
    if Zf[2] <= 5 and i < 100: #100 pt force averaging
        err.append(force)
        Rupture_ext = trap_ext
        Final_HBcount = Zf[2] #Zf[2] is HB count 
        i = i+1


y_peak = pd.Series(y).rolling(window=m_avg).mean().iloc[m_avg-1:].values # EMA calculation
y_pk = y_peak.tolist()
for loop in range(m_avg-1):
    y_pk.insert(0, 0)
    
break_force = np.mean(err)
error = np.std(err, ddof=1) / np.sqrt(np.size(err))
print (break_force, '+/-', error, 'pt-average-', i, 'Final hydrogenbond # ', Final_HBcount)

Loadrate = (ext_rate*l_units)/t_units
Loadrate = format(Loadrate, ".2e") + ' nm/s'
Force_label = format(break_force, ".2f") + ' pN'

#graph plotting 
ax.set(title='',ylabel='hydrogen bonds',xlabel='force (pN)')
ax.set(xlim=[00,120],ylim=[0,25])
plt.plot(x,y, color='#90BF2A',linewidth=0.45)
plt.plot(x,y_pk, color='#144C59',linewidth=0.075, alpha = 1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig('HB vs force for ' + str(i) +'_pt-avg_' + txt.rstrip('.dat') + '.png', transparent=False)
plt.show()
data.close()
