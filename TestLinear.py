import numpy as np
import matplotlib.pyplot as plt
import LossFunctions as LF
from scipy import interpolate
from scipy import stats


###################################
#  Options for Loss Calculations  #
###################################

plotLoss = True
saveLoss = "Q_cGTY.dat"
plotConvergence = False
plotDelta = False
plotV = False
plotrms = False
plotGamma = False
plotTau = False
plotTauRatio = True
plotMod = False

constG = True
constT = True
constY = True

##################################
#Parameters for Loss Calculations#
##################################
maxV = 0.1
maxT = 27
minT = 17
#vol = 21. * 21. * 21.
vol = 20.94*20.94*20.94
nTLS = 6.5
freq = 1e3

###################################
#  Staging for Loss Calculations  #
###################################
tlsdata = np.loadtxt("pp-tls.tot")
#LF.RemoveOutliers(tlsdata)
print tlsdata.shape
tlsdata = LF.RemoveDuplicates(tlsdata)
nPts = tlsdata.shape[0]
print tlsdata.shape
tMatrix = np.linspace(minT,maxT, num = 500)
print tMatrix
omega = 2 * np.pi * freq
density = nTLS / vol

#Import data
tlsdata[:,5] = 1.0 / tlsdata[:,5] * 0.0101804979 * 1e-12
tlsdata[:,6] = 1.0 / tlsdata[:,6] * 0.0101804979 * 1e-12
if constG:
    tlsdata[:,7] = np.mean(tlsdata[:, 7])
    #tlsdata[:, 7] = 1.7
if constY:
    tlsdata[:, 10] = np.mean(tlsdata[:,10])
if constT:
    tlsdata[:,5] = np.mean(tlsdata[:,5])
    tlsdata[:,6] = np.mean(tlsdata[:,6])
tlsdata[:,10] = tlsdata[:,10] * 6.3227e-7

print np.mean(tlsdata[:, 7]), np.mean(tlsdata[:, 8])

omegas = np.linspace(500,50000)
#omegas = [1118, 3377, 7307, 10077, 12152, 16880]

logOmega = []
inverseT = []
losses = []
for omega in omegas:
    omega = omega * 2 * np.pi
    loss = LF.calculateLossFunction(tlsdata[:], density, omega, tMatrix)
    losses.append(loss)
    logOmega.append(np.log(omega))
    inverseT.append(1.0 / tMatrix[np.argmax(loss)])
    print np.log(omega), tMatrix[np.argmax(loss)], np.amax(loss)

for i in xrange(0, len(omegas)):
    plt.plot(tMatrix, losses[i], label = str(omegas[i]) + " Hz")
    plt.plot(tMatrix[np.argmax(losses[i])], np.amax(losses[i]), 'go')
plt.legend()
plt.show()

m, b, R, Pval, stderr = stats.linregress(inverseT, y = logOmega)
X = np.linspace(inverseT[0], inverseT[-1])
Y = m * X[:] + b
plt.scatter(inverseT, logOmega)
plt.plot(X, Y)
#plt.show()
plt.savefig("Linear_Varying.png")

print m, b, R
kb = 8.617342E-5
print np.exp(-b), -m * kb
#print "Tau = " + '{6.4f}'.format(np.exp(-b)) + ", V = " + '{6.4f}'.format(-m * kb)
