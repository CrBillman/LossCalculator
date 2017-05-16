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

nArgs = len(sys.argv) 
if nArgs < 1:
    print "You need to specify a filename!"
    sys.exit()
fns = []
for i in xrange(1, nArgs):
    print i
    fns.append(sys.argv[i])
    print fns
for fn in fns:
    data = np.loadtxt(fn)
    plt.hist(data, bins = 50)
    plt.title("Distribution of Loss at 20 K")
    plt.xlabel(r'$Q^{-1}$', labelpad = -10)
    plt.ylabel("Number of Bootstrap Samples")
#plt.show()
plt.savefig("20kVariation.png")
