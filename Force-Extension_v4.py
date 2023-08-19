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
from scipy.signal import find_peaks


#Input file wtih time, energy, HBs, strand length and x, y, z co-ords for 2 the two traps, i.e. Zf[4:7], Zf[7:10]
folder = '/home/arventh/oxDNAdata/....' #input directory
txt = 'Observables.dat' #input file
os.chdir(folder)
data = open(txt, 'r')
i, x, y, y_peak = 0,[],[],[]

dt = 0.005 #simulation time units (Output should have time in MD units where MD units = steps * dt)
F_units = 48.63 #units in pN 
l_units = 0.8518 #units in nm 
t_units = 3.03e-12 #seconds
k1,k2 = 0.2,0.2 #simulation units

m_avg = 200 #exponential moving average smoothing of n-points
ext_rate = 0.2e-7 #rate of increase length (simulation units) per simulation step
F_dir = np.array([0,0,1]) #direction of trap movement
t1_init = np.array([0, 0, 0]) #moving trap
t2_init = np.array([0, 0, 0]) #fixed trap

UnitF_dir = F_dir/np.linalg.norm(F_dir) #unit vector
init_dist = np.dot((t1_init - t2_init), UnitF_dir) #initial distance between two traps along force dir
keff = ((k1*k2)/(k1+k2)) #1 unit of force constant (1 unit force/1 unit length) - 57.09 pN/nm

fig = plt.figure() 
ax = fig.add_subplot(1,1,1)

for line in data:
    Zf = np.array(line.split(), dtype=float)
    l_dna = Zf[3] #end to end strand length 
    trap1_current = t1_init + (Zf[0]/dt)*ext_rate*UnitF_dir #trap1 displacement update
    trap_ext1 = trap1_current - Zf[4:7] #vector subtraction of moving trap1 and attached nucleotide
    trap_ext2 = Zf[7:10]- t2_init #vector subtraction of fixed trap2 and attached nucleotide
    force = (np.dot((trap_ext1+trap_ext2),UnitF_dir))*keff*F_units #dot product of summation of 2 traps along force axis
    overall_ext = trap1_current - t2_init
    trap_ext = (np.dot(overall_ext,UnitF_dir)-init_dist)*l_units
    x.append(trap_ext) #this is extension of traps along force-axis not the DNA strand extension.
    y.append(force)
    if Zf[2] <= 4 and i == 0: #Zf[2] is HB count and this sets the range for peak picking
        Rupture_F = force
        Rupture_ext = trap_ext
        y_peak.append(force)
        i = 1

y_peak = pd.Series(y).rolling(window=m_avg).mean().iloc[m_avg-1:].values # EMA calculation
y_pk = y_peak.tolist()
for loop in range(m_avg-1):
    y_pk.insert(0, 0)

ax.set(ylabel='Force (pN)', xlabel='Extension (nm)')
ax.set(xlim=[0,18],ylim=[-5,60])
plt.plot(x,y, color='#D94862',linewidth=0.5)

Loadrate = (ext_rate*l_units)/t_units
Loadrate = format(Loadrate, ".2e") + ' nm/s'
Force_label = format(force, ".2f") + ' pN'
plt.plot(x,y_pk, color='#344973',linewidth=0.5, alpha = 1)

ax.legend([Loadrate, 'EMA ('+ str(m_avg)+ ')'], loc=0, frameon=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

peaks, ht = find_peaks(y_pk, height=10, prominence=3) #finds a list of peaks satisfying the criteria
for peakapeaka in peaks:
    if x[int(peakapeaka)] > Rupture_ext: #Range in which peaks are needed
        anno = str(round(y_pk[int(peakapeaka)], 2)) #+' pN' # peak annotation
        plt.plot(x[int(peakapeaka)],y_pk[int(peakapeaka)],"o", color='#000000', alpha=0.6, fillstyle = 'none')
        ax.annotate(anno,xy=(x[int(peakapeaka)], y_pk[int(peakapeaka)]), 
            xycoords='data', xytext=(0,0), #offset the rupture force annotation
            textcoords='offset points', horizontalalignment='left', verticalalignment='bottom')
 
for oolala in range(len(x)):
    csv = open("Graph_plot.csv", "a")
    csv.write(str(x[oolala]) + ',' + str(y[oolala]) + ',' + str(y_pk[oolala]) + '\n')

csv.close()
plt.show()
data.close()

