# -*- coding: utf-8 -*-
"""
plot and save the histgram of the X+ and X- for Gent and Rotterdam
Created on Wed Jan 28 11:43:14 2015

@author: weigang
"""

import numpy as np
import matplotlib.pylab as plt
import os

outputXGentGr = ['Ygr_'+str(n)+'.txt' for n in range(48)]
outputXGentSm = ['Ysm_'+str(n)+'.txt' for n in range(48)]
outputXRottGr = ['srcCateg1_Ygr.txt', 'srcCateg2_Ygr.txt', 'srcCateg3_Ygr.txt', 'srcCateg4_Ygr.txt']
outputXRottSm = ['srcCateg1_Ysm.txt', 'srcCateg2_Ysm.txt', 'srcCateg3_Ysm.txt', 'srcCateg4_Ysm.txt']

print os.path.exists("srcCateg1_Ygr.txt")

plt.figure(figsize=[12,5])
plt.subplot(1, 2, 1)
gr = []
# Gent input
for gent in outputXGentGr:
    print gent
    gr += list(np.loadtxt(gent))

# Rotterdam input
for rott in outputXRottGr:
    fo = open(rott, 'r')
    numlist = []
    print rott
    for n in range(300000): # load 300000 rows of the original file
        s = fo.readline()
        numlist.append(float(s))
    gr += numlist
    fo.close()
    
plt.hist(gr, np.linspace(0, 14, 42), normed=True, fill=True)
plt.ylabel('Percentage %')
plt.xlabel('Input variable')
plt.ylim([0,0.35])
plt.legend(['Xi+'])

plt.subplot(1, 2, 2)
sm = []
# Gent input
for gent2 in outputXGentSm:
    print gent2
    sm += list(np.loadtxt(gent2))

# Rotterdam input
for rott in outputXRottSm:
    fo = open(rott, 'r')
    numlist = []
    print rott
    for n in range(300000): # load 300000 rows of the original file
        s = fo.readline()
        numlist.append(float(s))
    sm += numlist
    fo.close()
    
plt.hist(sm, np.linspace(0, 14, 42), normed=True, fill=True, color='r')
plt.ylabel('Percentage %')
plt.xlabel('Input variable')
plt.legend(['Xi-'])
plt.ylim([0,0.35])
plt.show()