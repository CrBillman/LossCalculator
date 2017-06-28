import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import LossFunctions as LF
import Bootstrap as BS
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

plotLoss = True
bSamples = 3
percentiles = True
saveLoss = "Q_Percentiles1000.dat"
saveDistros = True

BoltzmannCorrection = False
constG = False
constT = False
constY = False

##################################
#Parameters for Loss Calculations#
##################################
maxV = 0.1
maxT = 300
minT = 1
nT = 100
vol = 20.94*20.94*20.94
nTLS = 7.0
freq = 1e3

###################################
#  Staging for Loss Calculations  #
###################################
omega = 2 * np.pi * freq
density = nTLS / vol
tMatrix = np.linspace(minT,maxT, num = nT)


if len(sys.argv) > 1:
    suff = "_" + sys.argv[1]
else:
    suff = ""


tlsdata = np.loadtxt("Data/pp-tls.tot" + suff)
tlsdata, nPts = LF.CleanData(tlsdata, constG = constG, constY = constY, constT = constT)

original = LF.calculateLossFunction(tlsdata, density, omega, tMatrix, BoltzmannCorrection = BoltzmannCorrection)

shuffleIndices = [0,3,4,5,6,7,10]
tlskde = BS.tlsKDE()
tlskde.fit(tlsdata, shuffleIndices, bandwidth = 0.186379686016)

losses = []
newtlsdata = np.copy(tlsdata)
for bi in xrange(0, bSamples):
    print "Bootstrapping: working on the {:d} step of {:d}.".format(bi + 1,bSamples)
    newtlsdata[:, shuffleIndices] = tlskde.sample(nPts)
    losses.append(LF.calculateLossFunction(newtlsdata, density, omega, tMatrix, BoltzmannCorrection = BoltzmannCorrection))

if saveDistros:
    BS.SaveDistros(tMatrix, original, losses)

lLoss, uLoss = BS.CalcConfidenceInterval(original, losses)
if plotLoss:
    LF.plotLoss(tMatrix, original, lLoss, uLoss)
if saveLoss:
    LF.saveLoss(tMatrix, original, lLoss, uLoss, saveLoss)
