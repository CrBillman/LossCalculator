import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn import neighbors
from sklearn.model_selection import GridSearchCV
import os

def trans(data, scale = False):
    data[:, 3] = np.log10(data[:, 3])
    data[:, 4] = np.log10(data[:, 4])
    if scale:
        ss = StandardScaler()
        data = ss.fit_transform(data)
        return data, ss
    else:
        return data

def detrans(data, ss = None):
    if ss:
        data = ss.inverse_transform(data)
    data[:, 3] = np.power(10, data[:, 3])
    data[:, 4] = np.power(10, data[:, 4])
    return data

def SaveDistros(tMatrix, original, losses):
    if not os.path.isdir("./Bootstrap_Distributions"):
        os.makedirs("./Bootstrap_Distributions")
    template = "Bootstrap_Distributions/LossDistro_{:d}K.dat"
    lMatrix = np.matrix(losses)
    for i in xrange(0, tMatrix.shape[0]):
        fn = template.format(int(tMatrix[i]))
        tempMatrix = np.copy(lMatrix[:, i])
        tempMatrix = np.transpose(np.insert(lMatrix[:, i], 0, original[i]))
        np.savetxt(fn, tempMatrix)


def CalcConfidenceInterval(original, losses):
    lMatrix = np.matrix(losses)
    nTemp = lMatrix.shape[1]
    nLosses = len(losses)
    tempMatrix = np.empty([nLosses])
    lower = np.empty([nTemp])
    upper = np.empty([nTemp])
    ave = np.empty([nTemp])
    for i in xrange(0, nTemp):
        for j in xrange(0, nLosses):
            tempMatrix[j] = lMatrix[j, i]
        lower[i] = 2 * original[i] - np.percentile(tempMatrix, 97.5)
        upper[i] = 2 * original[i] - np.percentile(tempMatrix, 2.5)
    return lower, upper
        
def GridSearchKDE(data):
    params = {'bandwidth': np.logspace(-3, 3, 50)}
    grid = GridSearchCV(neighbors.KernelDensity(), params)
    grid.fit(data)

    print("best bandwidth: {0}".format(grid.best_estimator_.bandwidth))

    params = {'bandwidth': np.linspace(-0.5, 0.5, 50) * grid.best_estimator_.bandwidth}
    grid = GridSearchCV(neighbors.KernelDensity(), params)
    grid.fit(data)

    print("best bandwidth: {0}".format(grid.best_estimator_.bandwidth))
    return grid.best_estimator_.bandwidth

def CleanData(data, constG = False, constY = False, constT = False, BoltzmannCorrection = True):
    data = RemoveOutliers(data)
    nPts = data.shape[0]
    print "Using " + str(nPts) + " TLS data points."
    if not BoltzmannCorrection:
        data = RemoveDuplicates(data)

    #Convert data
    data[:,5] = 1.0 / data[:,5] * 0.0101804979 * 1e-12
    data[:,6] = 1.0 / data[:,6] * 0.0101804979 * 1e-12
    data[:,[9,10, 11]] = data[:,[9, 10, 11]] * 6.3227e-7
    if constG:
        data[:,7] = np.mean(data[:, 7])
    if constY:
        data[:, [9,10, 11]] = np.mean(data[:,[9,10, 11]])
    if constT:
        data[:,5] = np.mean(data[:,5])
        data[:,6] = np.mean(data[:,6])

    return data, nPts

class tlsKDE():
    kde = None
    scaler = None

    def __init__(self):
        return

    def fit(self, data, indices, bandwidth = None):
        fitData = data[:, indices]
        fitData, self.scaler = trans(fitData, scale = True)
        if not bandwidth:
            bandwidth = GridSearchKDE(fitData)
        self.kde = neighbors.KernelDensity(bandwidth = bandwidth)
        self.kde.fit(fitData)

    def sample(self, nPts):
        outData = self.kde.sample(n_samples = nPts)
        outData = detrans(outData, self.scaler)
        return outData
