import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate as intp

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

def GetSmooth(T, Tnew, Q):
    tck = intp.splrep(T, Q, s=0.5)
    return intp.splev(Tnew, tck, der=0)


lw = 3.0

title1 = r'$Q^{-1}$'
title2 = r'$Q^{-1}_{\gamma}$'
title3 = r'$Q^{-1}_{\tau, \gamma}$'
title4 = r'$Q^{-1}_{IBS}$'
title5 = r'$Q^{-1}_{Bulk}$'

plt.figure(figsize = (14,8), dpi = 150)

plt.figtext(0.1025, 0.93, "a)", size = 'xx-large', weight = 'bold')
plt.figtext(0.375, 0.93, "b)", size = 'xx-large', weight = 'bold')
plt.figtext(0.65, 0.93, "c)", size = 'xx-large', weight = 'bold')

cGTY = np.loadtxt("Q_cGTY.dat")
cT = np.loadtxt("Q_cT.dat")
cN = np.loadtxt("Q_cN.dat")
Exp = np.loadtxt("Exp.dat")
Bulk = np.loadtxt("Bulk.dat")
Exp[:,1] = Exp[:,1] * 1e3
Bulk[:,1] = Bulk[:,1] * 0.1
#plot "Q_cGTY.dat" w l lw 2 lc rgb "blue" smooth csplines, \
#         "Q_cT.dat" w l lw 2 lc rgb "red" smooth csplines, \
#          "Q_cN.dat" w l lw 2 lc rgb "purple" smooth csplines, \
#           "Exp.dat" u 1:(1e3 * $2) w l lw 2 lc rgb "black", \
#            "Bulk.dat" u 1:(0.1 * $2) w l lt 10 lw 2 lc rgb "black"
plt.subplot(131)
plt.plot(cGTY[:,0], cGTY[:,1], label =  title1, lw = lw)
plt.plot(cT[:,0], cT[:,1], label = title2, lw = lw)
plt.plot(cN[:,0], cN[:,1], label = title3, lw = lw)
plt.plot(Exp[:,0], Exp[:,1], label = title4, lw = lw, color = 'k')
plt.plot(Bulk[:,0], Bulk[:,1], label = title5, lw = lw, color = 'k', ls = '--')
plt.xlabel("Temperature (K)", fontsize = 'x-large')
plt.ylabel(title1, fontsize = 'large', weight = 'bold', x = 1.05)
plt.gca().yaxis.set_label_coords(-0.11,0.5)
plt.legend(fontsize = 'large')

plt.subplot(132)
plt.plot(cGTY[:,0], cGTY[:,1], label =  title1, lw = lw)
plt.plot(cT[:,0], cT[:,1], label = title2, lw = lw)
plt.plot(cN[:,0], cN[:,1], label = title3, lw = lw)
plt.plot(Exp[:,0], Exp[:,1], label = title4, lw = lw, color = 'k')
plt.plot(Bulk[:,0], Bulk[:,1], label = title5, lw = lw, color = 'k', ls = '--')
plt.xlabel("Temperature (K)", fontsize = 'x-large')
plt.ylabel(title1, fontsize = 'large', weight = 'bold')
plt.gca().yaxis.set_label_coords(-0.11,0.5)
plt.xlim([5, 100])
plt.legend(fontsize = 'large')

Tnew = np.linspace(5, 100, num = 300)

plt.subplot(133)
plt.plot(Tnew, GetSmooth(cGTY[:,0], Tnew, cGTY[:,1]), label =  title1, lw = lw)
plt.plot(Tnew, GetSmooth(cT[:,0], Tnew, cT[:,1]), label = title2, lw = lw)
plt.plot(Tnew, GetSmooth(cN[:,0], Tnew, cN[:,1]), label = title3, lw = lw)
#plt.plot(Tnew, GetSmooth(Exp[:,0], Tnew, Exp[:,1]), label = title4, lw = lw, color = 'k')
#plt.plot(Tnew, GetSmooth(Bulk[:,0], Tnew, Bulk[:,1]), label = title5, lw = lw, color = 'k', ls = '--')
plt.plot(Exp[:,0], Exp[:,1], label = title4, lw = lw, color = 'k')
plt.plot(Bulk[:,0], Bulk[:,1], label = title5, lw = lw, color = 'k', ls = '--')
plt.xlabel("Temperature (K)", fontsize = 'x-large')
plt.ylabel(title1, fontsize = 'large', weight = 'bold', x = 1.05)
plt.gca().yaxis.set_label_coords(-0.11,0.5)
plt.xlim([5, 100])
plt.legend(fontsize = 'large')
#plt.show()
plt.savefig("constComps.eps")
