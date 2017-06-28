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

data = np.loadtxt('Q_Percentiles100.dat')

wavens = np.array([473, 633, 974, 1300, 1544]) * 1e-9
freqs = 2.99792458e8 / wavens
alphas = np.array([2.98, 2.935, 2.9, 2.878, 2.865])

data[:, [1,2]] = data[:, [1,2]] * 1e24

print freqs, alphas

plt.figure(figsize = (16, 10), dpi = 300)

plt.subplot(121)
plt.plot(data[:, 0] / (2 * np.pi), data[:, 1], label = "Calculated")
#plt.scatter(freqs, alphas, s = 50, label = "Experimental Background", color = 'orange', marker = '^', zorder = 3)
#plt.axhline(2.91, color = 'black', ls = '--', label = "Estimated Background")
#plt.axhline(3.17, color = 'gray', ls = '-.', label = "Calculated Background")
plt.title("Real Part of the Polarizability")
plt.ylabel(r'$\alpha^{\prime} (cm^3) \times 10^{-24}$')
plt.xlabel("Frequency (Hz)")
plt.xscale('log')
plt.xlim(min(data[:, 0] / (2 * np.pi)), max(data[:, 0] / (2 * np.pi)))
#plt.legend()
#plt.yscale('log')

plt.subplot(122)
plt.plot(data[:, 0] / (2 * np.pi), data[:, 2])
plt.title("Imaginary Part of the Polarizability")
plt.ylabel(r'$\alpha^{\prime\prime} (cm^3) \times 10^{-24}$')
plt.xlabel("Frequency (Hz)")
plt.xscale('log')
#plt.yscale('log')
#plt.tight_layout()
#plt.show()
plt.savefig('Polarizability.png')
