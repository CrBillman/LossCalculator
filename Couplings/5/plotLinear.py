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
    coups = []
    cCoup = []
    strains = []
    cStrain = []
    while True:
        line = f.readline()
        if not line:
            break
        sLine = line.split()
        if len(sLine) == 0:
            coups.append(cCoup)
            strains.append(cStrain)
            cCoup = []
            cStrain = []
        else:
            cCoup.append(float(sLine[2]) - float(sLine[1]))
            cStrain.append(float(sLine[0]))
    return np.array(coups), np.array(strains[0])

def GetSlopes(pref):
    f = open(pref + "cc_ec_final.dump", 'r')
    f.readline()
    return np.array(map(float, f.readline().split()[:6]))

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
        plt.scatter(strain, coups[i, :], label = labels[i], color = colors[i])
        Y = 2*np.sign(coups[i, -1] - coups[i, 0]) * slopes[i] * strain[:] + coups[i, 3]
        plt.plot(strain, Y, color = colors[i])
    plt.legend()
    #plt.show()

calcs = [3, 5, 7, 8, 9, 11, 12, 14, 15, 16, 18 ,19, 20, 21, 22, 23]
nCols = 3
nRows = np.ceil(float(len(calcs))/ 3.0)
for i in xrange(0, len(calcs)):
    pref = str(calcs[i]) + "/"
    plt.subplot(nRows, nCols, i + 1)
    couplings, strain = ImportData(pref + "cc_ec_raw_data.dump")
    slopes =  GetSlopes(pref)
    PlotLines(couplings, strain, slopes)
    plt.title("Calculation " + str(cals))
plt.show()
