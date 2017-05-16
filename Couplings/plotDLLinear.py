import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

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

def ImportData(fn):
    f = open(fn, 'r')
    f.readline()
    energies = []
    cE = []
    strains = []
    cStrain = []
    while True:
        line = f.readline()
        if not line:
            break
        sLine = line.split()
        if sLine[0] == "0.00500000000":
            cStrain.append(float(sLine[0]))
            cE.append(float(sLine[1]))
            energies.append(cE)
            strains.append(cStrain)
            cE = []
            cStrain = []
        else:
            cStrain.append(float(sLine[0]))
            cE.append(float(sLine[1]))
    return np.array(energies), np.array(strains[0])

def GetSlopes(i):
    f = open("result_slope", 'r')
    while True:
        line = f.readline()
        if not line:
            break
        sLine = line.split()
        print sLine[0], i
        if sLine[0] == str(i):
            return np.array(map(float, sLine[1:7]))
    return np.array(xrange(0, 7))

def PlotLines(coups, strain, slopes):
    labels = {
            0 : r'$\gamma_{xx}$',
            1 : r'$\gamma_{yy}$',
            2 : r'$\gamma_{zz}$',
            3 : r'$\gamma_{xy}$',
            4 : r'$\gamma_{xz}$',
            5 : r'$\gamma_{yz}$'
            }
    colors = {
            0 : 'blue',
            1 : 'red',
            2 : 'orange',
            3 : 'purple',
            4 : 'green',
            5 : 'yellow'
            }
    for i in xrange(0, coups.shape[0]):
        Y = 2*np.sign(coups[i, -1] - coups[i, 0]) * slopes[i] * strain[:] + coups[i, 10]
        plt.plot(strain, Y, color = colors[i])
        plt.scatter(strain, coups[i, :], color = 'black', zorder = 3, s = 30)
        plt.scatter(strain, coups[i, :], label = labels[i], color = colors[i], zorder = 4)
    #plt.legend()
    #plt.show()

calcs = [3, 5, 7, 8, 9, 11, 12, 14, 15, 16, 18 ,19, 20, 21, 22, 23]
#calcs = [3]
nCols = 3
nRows = np.ceil(float(len(calcs))/ 3.0)
fig = plt.figure(figsize=(16, 16), dpi = 150)
for i in xrange(0, len(calcs)):
    pref = str(calcs[i]) + "/"
    plt.subplot(nRows, nCols, i + 1)
    E1, strain = ImportData(pref + "min1_cc/eng.dat")
    E2, strain = ImportData(pref + "min2_cc/eng.dat")
    couplings = E2 - E1
    slopes =  GetSlopes(calcs[i])
    print couplings
    print slopes
    PlotLines(couplings, strain, slopes)
    plt.title("Calculation " + str(calcs[i]))
    if i == len(calcs) - 1:
        plt.legend(bbox_to_anchor=(3.00, 0.85), ncol = 3)
plt.tight_layout()
#plt.show()
plt.savefig("DLCCs.png")
