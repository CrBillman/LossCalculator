import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import LossFunctions as LF
from scipy import interpolate
from sklearn import neighbors
import sys

mpl.rcParams['lines.linewidth'] = 4
mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['axes.labelsize'] = 'xx-large'
mpl.rcParams['legend.fontsize'] = 'xx-large'
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'
mpl.rcParams['xtick.major.width'] = 2.5
mpl.rcParams['ytick.major.width'] = 2.5
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['axes.titlesize'] = 'xx-large'
mpl.rcParams['axes.titleweight'] = 'bold'


###################################
#  Options for Loss Calculations  #
###################################

plotLoss = False
bSamples = 100
percentiles = True
saveLoss = "Q_Percentiles100.dat"
plotConvergence = False
plotDelta = False
plotV = False
plotrms = False
plotGamma = False
plotTau = False
plotTauRatio = True
plotMod = False

constG = False
constT = False
constY = False

##################################
#Parameters for Loss Calculations#
##################################
maxV = 0.1
maxT = 300
minT = 1
#vol = 21. * 21. * 21.
vol = 20.94*20.94*20.94
nTLS = 7.0
freq = 1e3

###################################
#  Staging for Loss Calculations  #
###################################
if len(sys.argv) > 1:
    suff = "_" + sys.argv[1]
else:
    suff = ""

tlsdata = np.loadtxt("pp-tls.tot" + suff)
tlsdata = LF.RemoveOutliers(tlsdata)
print tlsdata.shape
#tlsdata = LF.RemoveDuplicates(tlsdata)
nPts = tlsdata.shape[0]
tMatrix = np.linspace(minT,maxT, num = 100)
omega = 2 * np.pi * freq
density = nTLS / vol

print np.mean(tlsdata[:, 7])

#Import data
tlsdata[:,5] = 1.0 / tlsdata[:,5] * 0.0101804979 * 1e-12
tlsdata[:,6] = 1.0 / tlsdata[:,6] * 0.0101804979 * 1e-12
if constG:
    tlsdata[:,7] = np.mean(tlsdata[:, 7])
    #tlsdata[:, 7] = 1.7
if constY:
    tlsdata[:, [9,10, 11]] = np.mean(tlsdata[:,[9,10, 11]])
if constT:
    tlsdata[:,5] = np.mean(tlsdata[:,5])
    tlsdata[:,6] = np.mean(tlsdata[:,6])
tlsdata[:,[9,10, 11]] = tlsdata[:,[9, 10, 11]] * 6.3227e-7

original = LF.calculateLossFunction(tlsdata, density, omega, tMatrix, BoltzmannCorrection = True)

tlsProps = tlsdata[:,[0,1,5,6,7,10]]
print tlsProps[2, :]
tlsProps = LF.trans(tlsProps)
print tlsProps[2, :]
#kde = neighbors.KernelDensity(bandwidth = 0.008)
kde = neighbors.KernelDensity(bandwidth = 0.0001)
kde.fit(tlsProps)
losses = []
plt.hist(tlsProps[:,0], label = "Old Props", bins = 50)
newtlsdata = np.copy(tlsdata)
plt.hist(tlsdata[:,7], bins = 50, label = "Original")
for bi in xrange(0, bSamples):
    newProps = kde.sample(n_samples = nPts)
    newtlsdata[:,[0,1,5,6,7,10]] = tlsProps[:]
    newProps = LF.detrans(newProps)
    #newProps = LF.detrans(tlsProps)
    print newProps[2, :]
    newtlsdata[:,[0,1,5,6,7,10]] = newProps[:]
    print "Working on the " + str(bi + 1) + " step of " + str(bSamples) + "."
    losses.append(LF.calculateLossFunction(newtlsdata, density, omega, tMatrix, BoltzmannCorrection = True))
    #plt.hist(newtlsdata[:,7], bins = 50, label = "Generated with KDE" + str(bi))
    #plt.hist(newtlsdata[:,5], bins = 50, label = "Generated with KDE" + str(bi))
    #plt.plot(tMatrix, losses[bi], label = "New Loss" + str(bi))
#plt.legend()
#plt.ylim([0,1])
#plt.xlim([1e-15, 1e-14])
lLoss, avLoss, uLoss = LF.CalcConfidenceInterval(original, losses, saveDistros = True, percentiles = percentiles)
if plotLoss:
    plt.plot(tMatrix, avLoss, color = 'blue', label = "Calculated Loss")
    plt.plot(tMatrix, lLoss, color = 'orange', label = "95% Confidence Interval", lw = 2, ls = '--')
    plt.plot(tMatrix, uLoss, color = 'orange', lw = 2, ls = '--')
    plt.legend()
    plt.ylim([0.0, 1.0])
    plt.show()
if saveLoss:
    prntMat = np.empty([tMatrix.shape[0], 4])
    prntMat[:, 0] = tMatrix[:]
    prntMat[:, 1] = avLoss[:]
    prntMat[:, 2] = lLoss[:]
    prntMat[:, 3] = uLoss[:]
    np.savetxt(saveLoss, prntMat)
