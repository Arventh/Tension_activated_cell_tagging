#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 15:23:29 2022

@author: arventh
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

folder = '/oxDNAdata/Rong/24mer_peeling' #raw_input directory
txt = 'Graph_plot.csv'
os.chdir(folder)
data = open(txt, 'r')

x1, y1, x2, y2, x3, y3, x4, y4 = [],[],[],[],[],[],[],[]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for line in data:
    Zf = line.split(',')
    x1.append(Zf[0])
    y1.append(Zf[1])
    x2.append(Zf[2]) 
    y2.append(Zf[3])
    x3.append(Zf[4])
    y3.append(Zf[5])
    x4.append(Zf[6])
    y4.append(Zf[7])

ax.set(title='',ylabel='hydrogen bonds',xlabel='force (pN)')
ax.set(xlim=[00,120],ylim=[0,24])
plt.plot(x1,y1, color='#90BF2A',linewidth=0.45)
plt.plot(x2,y2, color='#144C59',linewidth=0.075, alpha = 1)
plt.plot(x3,y3, color='#90BF2A',linewidth=0.45)
plt.plot(x4,y4, color='#144C59',linewidth=0.075, alpha = 1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig('HB vs force for ' + txt.rstrip('.dat') + '.png', transparent=False)
plt.show()
data.close()
