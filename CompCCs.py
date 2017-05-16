import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import LossFunctions as LF
from scipy import interpolate
from scipy import stats
import sys

mpl.rcParams['lines.linewidth'] = 4
mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['axes.labelsize'] = 'x-large'
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
saveLoss = "Q_cGTY.dat"

##################################
#Parameters for Loss Calculations#
##################################

###################################
#  Staging for Loss Calculations  #
###################################
if len(sys.argv) > 1:
    suff = "_" + sys.argv[1]
else:
    suff = ""

tlsdata = np.loadtxt("pp-tls.tot" + suff)
#tlsdata = LF.RemoveDuplicates(tlsdata)
tlsdata = LF.RemoveOutliers(tlsdata)
print tlsdata.shape
nPts = tlsdata.shape[0]

#Import data
tlsdata[:,5] = 1.0 / tlsdata[:,5] * 0.0101804979 * 1e-12
tlsdata[:,6] = 1.0 / tlsdata[:,6] * 0.0101804979 * 1e-12
tlsdata[:, [7,8] ] = tlsdata[:, [7,8] ] * 2.0
tlsdata[:,10] = tlsdata[:,10] * 6.3227e-7

fig = plt.figure(figsize = (16, 12), dpi = 150)
#fig = plt.figure(dpi = 150)
plt.subplot(1, 2, 1)
plt.title("Comparison of Coupling Constant Distributions")
plt.hist(tlsdata[:,7], bins = 100, label = "Sim. Long.")
plt.axvline(1.04, color = 'blue', lw = 6, ls = '--', label = "Exp. Long.")
plt.hist(tlsdata[:,8], bins = 50, label = "Sim. Trans.")
plt.axvline(0.65, color = 'green', lw = 6, ls = '--', label = "Exp. Trans.")
plt.legend()
plt.xlim([0,2])
plt.xlabel("Coupling Constant (eV)")
plt.ylabel("Number of TLS's")

print np.mean(tlsdata[:, 7]), 1.04

plt.subplot(1, 2, 2)
gRatio = np.divide(tlsdata[:,7], tlsdata[:,8])
plt.title("Ratio of Coupling Constants")
plt.hist(gRatio, bins = 150, label = "Simulated")
plt.xlabel("Ratio of " + r'$\gamma_L :\gamma_T$')
plt.ylabel("Number of TLS's")
plt.axvline(np.mean(gRatio), label = "Mean", lw = 6, color = 'orange')
plt.axvline(np.percentile(gRatio, 50), label = "Median", lw = 6, color = 'red')
plt.axvline(np.sqrt(2.5), label = "Exp.", lw = 6, color = 'black')
plt.xlim([0, 20])
plt.legend()
#plt.show()
plt.tight_layout()
plt.savefig("CompCCs.png")

print np.sqrt(2.5), np.percentile(gRatio, 50), np.mean(gRatio)
