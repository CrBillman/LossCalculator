import numpy as np
import matplotlib.pyplot as plt
import LossFunctions as LF

def checkFile(fname):
	nCols = 12
	stream = open(fname, "r")
	fixed = open("tls.checked", "wb")
	i = 1
	while True:
		line = stream.readline()
		if(not line):
			break
		sLine = line.split()
		if(len(sLine) != nCols):
			print i
			print sLine
		else:
			fixed.write(line)
		i = i + 1
	return


BoltzmannFactor = False

checkFile("tls.tot")
data = np.loadtxt("pp-tls.tot")
print "About to"
data = LF.RemoveDuplicates(data)
#data = LF.RemoveOutliers(data)
print max(data[:,7])
print "Did it"
print data.shape

data[:,5] = 1.0 / data[:,5] * 0.0101804979
data[:,10] = data[:,10] * 6.3227e-7

AvgTime = 0.0
AvgLCoupling = 0.0
AvgTCoupling = 0.0
for i in xrange(0, data.shape[0]):
	AvgTime = AvgTime + data[i,5] * np.exp(data[i,1] / 0.086173)
	AvgLCoupling = AvgLCoupling + data[i,7] * np.exp(data[i,1] / 0.086173)
	AvgTCoupling = AvgTCoupling + data[i,8] * np.exp(data[i,1] / 0.086173)

print "Average Relaxation Time (ps): ", AvgTime / data.shape[0], np.mean(data[:,5])
print "Average Longitudinal Coupling Constant (eV/A): ", AvgLCoupling / data.shape[0], np.mean(data[:,7])
print "Average Transverse Coupling Constant (eV/A): ", AvgTCoupling / data.shape[0], np.mean(data[:,8])
print "Average Bulk Modulus (eV/A^3): ", np.mean(data[:,10]), np.std(data[:,10])
