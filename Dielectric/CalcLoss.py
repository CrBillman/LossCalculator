import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import LossFunctions as LF
import DielectricLoss as DL
from scipy import interpolate
from sklearn import neighbors
import sys

kb = 8.617342E-5

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
BoltzmannCorrection = False

bootStrap = False
constT = False
dipole = 3
temp = 300

dipole = dipole * 3.33564e-30

##################################
#Parameters for Loss Calculations#
##################################
maxV = 0.1
maxF = 1e17
minF = 1e11
nF = 1e2
#vol = 21. * 21. * 21.
#vol = 20.94*20.94*20.94
vol = 1.0
nTLS = 1.0
#freq = 1e3

#vol = vol * 1e-30

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
fMatrix = np.logspace(np.log10(minF), np.log10(maxF), num = nF)
fMatrix = fMatrix[:] * 2 * np.pi
density = nTLS / vol

print np.mean(tlsdata[:, 7])

tlsdata[:,5] = 1.0 / tlsdata[:,5] * 0.0101804979 * 1e-12
tlsdata[:,6] = 1.0 / tlsdata[:,6] * 0.0101804979 * 1e-12
tlsdata[:,7] = dipole
tlsdata[:,8] = dipole

originalLoss = DL.calculateSusceptibility(tlsdata, density, temp, fMatrix, BoltzmannCorrection = BoltzmannCorrection)
originalLoss = originalLoss[:] * 8.988e15

if saveLoss:
    prntMat = np.empty([fMatrix.shape[0], 3])
    prntMat[:, 0] = fMatrix[:]
    prntMat[:, 1] = np.real(originalLoss[:])
    prntMat[:, 2] = np.imag(originalLoss[:])
    np.savetxt(saveLoss, prntMat)
