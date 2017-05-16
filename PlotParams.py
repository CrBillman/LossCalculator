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
mpl.rcParams['axes.labelsize'] = 'x-large'
mpl.rcParams['legend.fontsize'] = 'xx-large'
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'
mpl.rcParams['xtick.major.width'] = 2.5
mpl.rcParams['ytick.major.width'] = 2.5
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['axes.titlesize'] = 'xx-large'
mpl.rcParams['axes.titleweight'] = 'bold'

def AdjustXAxis(data):
    ll = np.percentile(data, 2)
    ul = np.percentile(data, 98)
    plt.xlim([ll, ul])
    plt.gca().xaxis.grid(True)
    return


###################################
#  Options for Loss Calculations  #
###################################
plotDupes = False
mod = 11
##################################
#Parameters for Loss Calculations#
##################################
maxV = 0.1
maxT = 300
minT = 1
#vol = 21. * 21. * 21.
vol = 20.94*20.94*20.94
nTLS = 6.5
freq = 1e3

###################################
#  Staging for Loss Calculations  #
###################################
if len(sys.argv) > 1:
    suff = "_" + sys.argv[1]
else:
    suff = ""

tlsdata = np.loadtxt("pp-tls.tot" + suff)
#LF.RemoveOutliers(tlsdata)
nPts = tlsdata.shape[0]

#Import data
tlsdata[:,5] = 1.0 / tlsdata[:,5] * 0.0101804979
tlsdata[:,6] = 1.0 / tlsdata[:,6] * 0.0101804979
tlsdata[:,9] = tlsdata[:,9] * 6.3227e-7
tlsdata[:,10] = tlsdata[:,10] * 6.3227e-7
tlsdata[:,11] = tlsdata[:,11] * 6.3227e-7

NoDupes = LF.RemoveDuplicates(tlsdata)

#q = data[np.percentile(data[:, 7] , 99), 7]
grt = np.divide(tlsdata[:,7], tlsdata[:,8])
grtSq = np.multiply(grt, grt)
print "Average Ratio of Coupling Constants Squared is: " + str(np.mean(grtSq))

print np.mean(tlsdata[:, 7]), np.mean(tlsdata[:, 8])

gcf = plt.figure(figsize = (16,12), dpi = 200)
plt.figtext(0.0525, 0.91, "A", size = 'xx-large', weight = 'bold')
plt.figtext(0.0525, 0.64, "B", size = 'xx-large', weight = 'bold')
plt.figtext(0.0525, 0.37, "C", size = 'xx-large', weight = 'bold')
plt.figtext(0.4825, 0.91, "D", size = 'xx-large', weight = 'bold')
plt.figtext(0.4825, 0.64, "E", size = 'xx-large', weight = 'bold')
plt.figtext(0.4825, 0.37, "F", size = 'xx-large', weight = 'bold')
plt.subplot2grid((12,1),(0,0), rowspan = 3)

plt.subplot2grid((9,11), (0, 0), rowspan = 2, colspan = 5)
plt.hist(tlsdata[:,1], bins = 150, label = "All points", normed = True)
if plotDupes:
    plt.hist(NoDupes[:,1], bins = 150, label = "After Duplicate Removal", normed = True)
    plt.legend(fontsize = 'medium')
AdjustXAxis(tlsdata[:,1])
plt.xlabel("Barrier Height (eV)")
plt.ylabel("Density (Arb. units)")
plt.title("Barrier Height")

plt.subplot2grid((9,11), (0, 6), rowspan = 2, colspan = 5)
plt.hist(tlsdata[:,0], bins = 200, label = "All points", normed = True)
if plotDupes:
    plt.hist(NoDupes[:,0], bins = 200, label = "After Duplicate Removal", normed = True)
    plt.legend(fontsize = 'medium')
AdjustXAxis(tlsdata[:,0])
#plt.xlim([0,0.06])
plt.xlabel("Asymmetry (eV)")
plt.ylabel("Density (Arb. units)")
plt.title("Asymmetry")

plt.subplot2grid((9,11), (3, 0), rowspan = 2, colspan = 5)
plt.hist(tlsdata[:,2], bins = 75, label = "All points", normed = True)
if plotDupes:
    plt.hist(NoDupes[:,2], bins = 75, label = "After Duplicate Removal", normed = True)
    plt.legend(fontsize = 'medium')
AdjustXAxis(tlsdata[:,2])
plt.xlabel("rms distance (Ang.)")
plt.ylabel("Density (Arb. units)")
plt.title("RMS Distance")

plt.subplot2grid((9,11), (3, 6), rowspan = 2, colspan = 5)
plt.hist(np.log10(tlsdata[:,5]), bins = 150, label = "All points", normed = True)
if plotDupes:
    plt.hist(np.log10(NoDupes[:,5]), bins = 150, label = "After Duplicate Removal", normed = True)
    plt.legend(fontsize = 'medium', loc = 2)
AdjustXAxis(np.log10(tlsdata[:,5]))
plt.xticks([-4, -3, -2, -1], [1e-4, 1e-3, 1e-2, 1e-1])
plt.xlabel("Relaxation Time (ps)")
plt.ylabel("Density (Arb. units)")
plt.title("Average Relaxation Time")

plt.subplot2grid((9,11), (6, 0), rowspan = 2, colspan = 5)
plt.hist(tlsdata[:,7], bins = 150, label = "All points", normed = True)
if plotDupes:
    plt.hist(NoDupes[:,7], bins = 150, label = "After Duplicate Removal", normed = True)
    plt.legend(fontsize = 'medium')
AdjustXAxis(tlsdata[:,7])
plt.xlabel("Coupling Constant (eV)")
plt.ylabel("Density (Arb. units)")
plt.title("Longitudinal Coupling Constant")

plt.subplot2grid((9,11), (6, 6), rowspan = 2, colspan = 5)
avg = np.mean(tlsdata[:,mod])
plt.hist(tlsdata[:,mod], bins = 150, label = "All points", normed = True)
if plotDupes:
    plt.hist(NoDupes[:,mod], bins = 150, label = "After Duplicate Removal", normed = True)
    plt.legend(fontsize = 'medium', loc = 2)
AdjustXAxis(tlsdata[:,mod])
plt.ylabel("Density (Arb. units)")
if mod == 9:
    modStr = "Bulk"
elif mod == 10:
    modStr = "Shear"
elif mod == 11:
    modStr = "Young's"
plt.xlabel(modStr + " Modulus (eV)")
plt.title(modStr + " Modulus")
#plt.show()
plt.savefig("TLS_Stats" + suff + ".eps")
print avg, np.std(tlsdata[:,mod]), max(tlsdata[:, mod]), min(tlsdata[:, mod])
