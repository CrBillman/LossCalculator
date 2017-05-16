#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 15:44:29 2017

@author: Chris
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import LossFunctionsp3 as LF
from scipy import stats
from scipy import interpolate

tlsdata = np.loadtxt("pp-tls.tot")
LF.RemoveOutliers(tlsdata)
tlsdata = LF.RemoveDuplicates(tlsdata)
tlsdata = LF.RemoveOutliers(tlsdata)
tlsdata[:,5] = 1.0 / tlsdata[:,5] * 0.0101804979# * 1e-12
tlsdata[:,6] = 1.0 / tlsdata[:,6] * 0.0101804979# * 1e-12
tlsdata[:,10] = tlsdata[:,10] * 6.3227e-7
       
nPts = tlsdata.shape[0]

tlsdata[:, [5,6]] = np.log10(tlsdata[:, [5,6]])
#sGV = np.sort(np.asarray(tlsdata)[:,[1,])
#i2p5 = int(0.025 * nPts)
#i97p5 = int(0.975 * nPts)
#ssV = sV[i2p5:i97p5]
#ssG = sG[i2p5:i97p5]
ssV = tlsdata[:,1]
ssD = tlsdata[:,0]
ssG = tlsdata[:, 7]
fig = plt.figure(figsize = (20,20), dpi = 150)
indices = [0, 1, 2, 5, 6, 7, 10]
#title = ["Asymmetry", "Barrier Height", "RMS Distance", "Tau 1", "Tau 2", 
#         "Coupling Constant", "Modulus"]
title = [r'$\Delta$', r'$V$', r'$d_{RMS}$', r'$\tau_1$', r'$\tau_2$',
         r'$\gamma_L$', r'$\epsilon$']
xshift = [0, 0, -0.015, 0, 0, 0, 0]
q2 = [85, 95, 95, 95, 95, 95, 95]
ppr = len(indices)
fx = 0.935
ix = 0.165
fy = 0.815
iy = 0.065
ly = fy - iy
dy = ly / float(ppr -1)
lx = fx - ix
dx = lx / float(ppr-1)

pw = 1
ph = 1
ci = pw + 1
ri = ph + 1

for i in range(0, ppr):
    for j in range(0, ppr):
        ii = indices[i]
        jj = indices[j]
        #plt.subplot2grid((ci*ppr+1, ri*ppr+1), (ci*i+1, ri*j+1), colspan = pw, rowspan = ph)
        plt.subplot2grid((ppr+1, ppr+1), (i+1, j+1), colspan = pw, rowspan = ph)
        #plt.hist2d(np.asarray(tlsdata)[:,i], np.asarray(tlsdata)[:,j], bins = 200, norm = LogNorm())
        plt.hist2d(np.asarray(tlsdata)[:,jj], np.asarray(tlsdata)[:,ii], bins = 100)
        x1 = np.percentile(tlsdata[:,jj], 5)
        x2 = np.percentile(tlsdata[:,jj], q2[j])
        y1 = np.percentile(tlsdata[:,ii], 5)
        y2 = np.percentile(tlsdata[:,ii], q2[i])
        plt.xlim([x1, x2])
        plt.ylim([y1, y2])
        if ii == 5 or ii == 6:
            plt.yticks([-4, -3, -2, -1], [r'$10^{-16}$', r'$10^{-15}$', r'$10^{-14}$', r'$10^{-13}$'])
        if jj == 5 or jj == 6:
            plt.xticks([-4, -3, -2, -1], [r'$10^{-16}$', r'$10^{-15}$', r'$10^{-14}$', r'$10^{-13}$'])
for i in range(0, ppr):
    fig.text(ix + dx * i + xshift[i], 0.89, title[i], weight = 'bold', fontsize = 40, zorder = 5)
    fig.text(0.06 + xshift[-1-i], iy + dy * i, title[-1-i], weight = 'bold', fontsize = 40, zorder = 5)
    #plt.text(0.5, 0.5, title[i])
plt.tight_layout()
plt.savefig("test.png")
#print(ssV)
#print(stats.pearsonr(ssD, ssG))
#print(stats.spearmanr(ssD, ssG))
#plt.hist2d(np.asarray(tlsdata)[:,0], np.asarray(tlsdata)[:,7], bins = 200, norm = LogNorm())
#plt.hist2d(np.asarray(tlsdata)[:,0], np.asarray(tlsdata)[:,1], bins = 200, norm = LogNorm())
#plt.colorbar()
#plt.plot([0, 0.1], [0, 0.05], color = "red", lw = 4)
#plt.xlim([0,0.1])
#plt.ylim([0,0.2])
#plt.show()