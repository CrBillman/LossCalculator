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

scale = 1.0
#scale = 0.15

cMatrix = np.loadtxt("Q_Percentiles100.dat")
cMatrix[:, [1, 2, 3]] = cMatrix[:, [1, 2, 3]] * scale
expMatrix = np.loadtxt("Exp.dat")
expMatrix[:, 1] = expMatrix[:, 1] * 1e3

plt.title("Scale of " + str(scale) + " Applied to Calculations")
plt.fill_between(cMatrix[:, 0], cMatrix[:, 2], cMatrix[:, 3], color = 'blue', label = "95% Confidence Interval")
plt.plot(cMatrix[:, 0], cMatrix[:, 1], color = 'orange', label = "Average")
plt.plot(expMatrix[:, 0], expMatrix[:, 1], label = "Experiment", color = 'black')
plt.legend()
plt.xlabel("Temperature (K)")
plt.ylabel(r'$Q^{-1}$')
plt.ylim([0,1.5])
plt.show()
