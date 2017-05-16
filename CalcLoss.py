import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import LossFunctions as LF
from scipy import interpolate
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
saveLoss = "Qnodupes.dat"
plotConvergence = True
plotDelta = False
plotV = False
plotrms = False
plotGamma = False
plotTau = False
plotTauRatio = False
plotMod = False

DL = False
BoltzmannCorrection = True

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
#nTLS = 6.5
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
#print "Code is removing conversions for a test!"
print tlsdata.shape
nPts = tlsdata.shape[0]
tMatrix = np.linspace(minT,maxT, num = 100)
omega = 2 * np.pi * freq
density = nTLS / vol

#Import data
if not DL:
    tlsdata[:,5] = 1.0 / tlsdata[:,5] * 0.0101804979 * 1e-12
    tlsdata[:,6] = 1.0 / tlsdata[:,6] * 0.0101804979 * 1e-12
    tlsdata[:,[9, 10, 11]] = tlsdata[:,[9, 10, 11]] * 6.3227e-7
if not BoltzmannCorrection:
    #tlsdata = LF.RemoveOutliers(tlsdata)
    tlsdata = LF.RemoveDuplicates(tlsdata)
if constG:
    tlsdata[:,7] = np.mean(tlsdata[:, 7])
    #tlsdata[:, 7] = 1.7
if constY:
    tlsdata[:, [9, 10, 11]] = np.mean(tlsdata[:,[9, 10, 11]])
if constT:
    tlsdata[:,5] = np.mean(tlsdata[:,5])
    tlsdata[:,6] = np.mean(tlsdata[:,6])

print "Average Young's modulus: ", np.mean(tlsdata[:, [11]])

if plotLoss:
    loss = LF.calculateLossFunction(tlsdata[:], density, omega, tMatrix, BoltzmannCorrection = BoltzmannCorrection)
    plt.plot(tMatrix, loss)
    plt.show()
    if saveLoss:
        prntMat = np.empty([loss.shape[0], 2])
        prntMat[:, 0] = tMatrix[:]
        prntMat[:, 1] = loss[:]
        np.savetxt(saveLoss, prntMat)

if plotConvergence:
    p95 = int(nPts * 0.95)
    p90 = int(nPts * 0.90)
    p85 = int(nPts * 0.85)
    if plotV:
        plt.hist(tlsdata[:,1], bins = 250, label = "All points")
        plt.hist(tlsdata[:p95,1], bins = 250, label = "95% of points")
        plt.hist(tlsdata[:p90,1], bins = 250, label = "90% of points")
        plt.hist(tlsdata[:p85,1], bins = 250, label = "85% of points")
        plt.legend()
        plt.xlim([0,0.15])
        plt.xlabel("Barrier Height (eV)")
        plt.ylabel("Number of TLS's")
        plt.title("g(V)")
        plt.show()

    if plotDelta:
        plt.hist(tlsdata[:,0], bins = 250, label = "All points")
        plt.hist(tlsdata[:p95,0], bins = 250, label = "95% of points")
        plt.hist(tlsdata[:p90,0], bins = 250, label = "90% of points")
        plt.hist(tlsdata[:p85,0], bins = 250, label = "85% of points")
        plt.legend()
        plt.xlim([0,0.06])
        plt.xlabel("Asymmetry (eV)")
        plt.ylabel("Number of TLS's")
        plt.title("f(delta)")
        plt.show()

    if plotrms:
        plt.hist(tlsdata[:,2], bins = 75, label = "All points")
        plt.hist(tlsdata[:p95,2], bins = 250, label = "95% of points")
        plt.hist(tlsdata[:p90,2], bins = 250, label = "90% of points")
        plt.hist(tlsdata[:p85,2], bins = 250, label = "85% of points")
        plt.legend()
        #plt.xlim([0,0.15])
        plt.xlabel("rms distance (Ang.)")
        plt.ylabel("Number of TLS's")
        plt.show()

    if plotTau:
        plt.hist(tlsdata[:,5], bins = 250, label = "All points")
        plt.hist(tlsdata[:p95,5], bins = 250, label = "95% of points")
        plt.hist(tlsdata[:p90,5], bins = 250, label = "90% of points")
        plt.hist(tlsdata[:p85,5], bins = 250, label = "85% of points")
        plt.legend()
        plt.xlim([0,0.51e-12])
        plt.xlabel("Relaxation Time (s)")
        plt.ylabel("Number of TLS's")
        plt.show()

    if plotTauRatio:
        error = (
            (0.05, 14.1124743036),
            (0.1, 12.2037046631),
            (0.125, 11.0974050379),
            (0.2, 8.11091903323),
            (0.5, 2.00415017792),
            (1.0, 0.000105844169534),
            (2.0, 2.08792689676),
            (5.0, 9.35874848105),
            (8.0, 13.585899382),
            (10.0, 15.4390564741),
            (20.0, 20.5263960811)
        )
        #errorArray = np.array(error)
        errorArray = np.loadtxt("errorArray.dat")
        errorArray = errorArray[errorArray[:,0].argsort()]
        errorArray[:, 0] = np.log10(errorArray[:,0])
        tck, u = interpolate.splprep(errorArray.transpose(), s = 1)
        unew = np.arange(0, 1.01, 0.01)
        out = interpolate.splev(unew, tck)

        tauRatio = tlsdata[:,5] / tlsdata[:,6]
        print tauRatio
        sumOutside = 0
        for i in xrange(0, tauRatio.shape[0]):
            if tauRatio[i] <= 0.1 or tauRatio[i] >= 10:
                sumOutside = sumOutside + 1
        print float(sumOutside) / float(tauRatio.shape[0])
        tauRatio = abs(np.log10(tauRatio[:]))
        plt.hist(tauRatio[:], bins = 50)
        #plt.hist(tauRatio[:p95], bins = 150, label = "95% of points")
        #plt.hist(tauRatio[:p90], bins = 150, label = "90% of points")
        #plt.hist(tauRatio[:p85], bins = 150, label = "85% of points")
        plt.xticks([1.0, 2.0, 3.0], [r'$10^1$', r'$10^2$', r'$10^3$'])
        ax1 = plt.gca()
        ax2 = ax1.twinx()
        #ax2.plot(out[0], out[1], color='orange')
        #ax2.plot(errorArray[:,0], errorArray[:, 1], 'go', label = "Calculated Error using Numerical Model")
        ax2.hist(tauRatio[:], bins=50, normed=1, histtype='step', cumulative=1, label='Cumulative', color = 'orange', lw = 4)
        ax1.set_xlabel("Ratio of Relaxation Times")
        ax1.set_ylabel("Density of TLS's (Arb. Units)")
        ax2.set_ylabel("Cumulative Percent")
        ax2.set_yticks([0.0, 0.25, 0.50, 0.75, 1.0])
        ax2.set_yticklabels([0, 25, 50, 75, 100])
        ax2.set_ylim([0.0, 1.05])
        #ax2.legend()
        #plt.legend()
        plt.xlim([min(tauRatio), max(tauRatio)])
        #plt.show()
        plt.savefig("HistRatioTaus.png", dpi = 150)

    if plotGamma:
        plt.hist(tlsdata[:,7], bins = 250, label = "All points")
        plt.hist(tlsdata[:p95,7], bins = 250, label = "95% of points")
        plt.hist(tlsdata[:p90,7], bins = 250, label = "90% of points")
        plt.hist(tlsdata[:p85,7], bins = 250, label = "85% of points")
        plt.legend()
        #plt.xlim([0,5])
        plt.xlabel("Coupling Constant (eV)")
        plt.ylabel("Number of TLS's")
        plt.show()

    if plotMod:
        exp = 0.225
        avg = np.mean(tlsdata[:,10])
        plt.hist(tlsdata[:,10], bins = 150, label = "Calculated from TLS's")
        plt.axvline(exp, color = 'red', alpha = 0.8, label = "Experimental Value", lw = 2)
        plt.axvline(avg, color = 'green', alpha = 0.8, label = "Average Calculated", lw = 2)
        print('{:25}'.format("Experiment (eV)") + '{:25}'.format("Average Calculated (eV)") + '{:25}'.format("Percent Error"))
        print('{:25}'.format('{:5.3f}'.format(exp)) + '{:25}'.format('{:5.3f}'.format(avg)) + '{:25}'.format('{:5.3f}'.format(abs(exp - avg) / exp * 100) + '%'))
        #plt.hist(tlsdata[:p95,10], bins = 250, label = "95% of points")
        #plt.hist(tlsdata[:p90,10], bins = 250, label = "90% of points")
        #plt.hist(tlsdata[:p85,10], bins = 250, label = "85% of points")
        plt.legend(loc = 2)
        #plt.xlim([0,5])
        plt.xlabel("Young's Modulus (eV)")
        plt.ylabel("Number of TLS's")
        plt.title("Young's Modulus Calculated for Silica Two-Level Systems")
        #plt.show()
