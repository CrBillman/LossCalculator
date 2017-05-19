import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn import neighbors
from sklearn.model_selection import GridSearchCV
import sys

def Tau(Tau1, Tau2, V, T, delta):
    numCut = 80
    kb = 8.617342E-5
    if (V+delta) / (kb * T)  > numCut or (V - delta) / (kb * T) > numCut:
        return np.exp(numCut)
    else:
        return 0.5 * (Tau1 * np.exp((V+delta)/(kb * T)) + Tau2 * np.exp((V-delta)/(kb * T)))

def resonanceFraction(Tau1, Tau2, V, T, delta, omega):
    tau = Tau(Tau1, Tau2, V, T, delta)
    if(np.isnan(tau)):
        return 0.0
    else:
        return tau * omega / (1 + omega **2 * tau **2)

def sechTerm(delta, delta0, T):
    numCut = 200
    kb = 8.617342E-5
    if np.sqrt(delta**2 + delta0**2)/(2 * kb * T) > numCut:
        return 0.0
    else:
        coshTerm = np.cosh(np.sqrt(delta**2 + delta0**2)/(2 * kb * T)) ** 2
    return 1.0/coshTerm

def dipoleMoment(delta, delta0, gamma):
    return delta**2 / (delta**2 + delta0**2) * gamma **2

def prefactor(Modulus, T):
    kb = 8.617342E-5
    return 1/(3.0 * kb * T * Modulus)

def calculateLossFunction(data, density, omega, tMatrix, BoltzmannCorrection = False):
    kb = 8.617342E-5
    qMatrix = np.zeros(tMatrix.shape)
    nTemp = tMatrix.shape[0]
    try:
        nPts = data.shape[0]
    except AttributeError:
        nPts = len(data)
    if BoltzmannCorrection:
        BCc = 0
        norm = 0.0
        for ip in xrange(0, len(data)):
            V1 = data[ip, 3]
            norm = norm + np.exp( (V1) / (kb * 1000))
        #norm = density / (norm)
        norm = density / float(nPts)
        print norm, density / float(nPts)
    else:
        BCc = 1
        norm = density / float(nPts)
    
    for it in xrange(0, nTemp):
        for ip in xrange(0,len(data)):
            delta = data[ip, 0]
            #V = data[ip, 1]
            V = 0.5 * (data[ip, 3] + data[ip, 4])
            V1 = data[ip, 3]
            tau1 = data[ip, 5]
            tau2 = data[ip, 6]
            Y = data[ip, 11]
            gamma = data[ip, 7]
            BC = max(BCc, np.exp( (V1) / (kb * 1000)) )
            qMatrix[it] = qMatrix[it] + BC * prefactor(Y,tMatrix[it]) * resonanceFraction(tau1, tau2, V, tMatrix[it], delta, omega) * sechTerm(delta, 1e-4, tMatrix[it]) * dipoleMoment(delta, 1e-4, gamma)
    qMatrix[:] = qMatrix[:] * norm * 1e3
    return qMatrix

def RemoveOutliers(data):
    qu = np.percentile(data[:, 7], 99)
    ql = np.percentile(data[:,7], 1)
    noOuts = []
    for i in xrange(0, data.shape[0]):
        gamma = data[i, 7]
        if gamma > ql and gamma < qu:
            noOuts.append(data[i, :])
    return np.matrix(noOuts)

def RemoveDuplicates(data):
    window = 200
    colsToComp = [0, 2, 3, 4]
    eComp = 0.005
    dComp = 0.001
    #colsToComp = [0, 1]
    cut = 0.00
    absCut = {0 : eComp, 
            2 : dComp, 
            3 : eComp,
            4 : eComp}
    nRows = data.shape[0]
    noDups = []
    for i in xrange(0, nRows):
        match = False
        #for j in xrange(i + 1, nRows):
        ul = min(i + window, nRows)
        for j in xrange(i + 1, ul):
            matchN = 1.0
            for k in colsToComp:
                if abs(data[i, k] - data[j, k]) <= absCut[k]:
                    matchN = matchN * 1.0
                else:
                    matchN = matchN * 0.0
            if matchN == 1.0:
                match = True
                break
        if not match:
            noDups.append(data[i,:])
            #print data[i, :]
    return np.matrix(noDups)

def trans(data, scale = False):
    #data[:, 2] = np.log10(data[:, 2])
    #data[:, 3] = np.log10(data[:, 3])
    data[:, 3] = np.log10(data[:, 3])
    data[:, 4] = np.log10(data[:, 4])
    if scale:
        ss = StandardScaler()
        data = ss.fit_transform(data)
        return data, ss
    else:
        return data

def detrans(data, ss = None):
    #data[:, 2] = np.power(10, data[:, 2])
    #data[:, 3] = np.power(10, data[:, 3])
    if ss:
        data = ss.inverse_transform(data)
    data[:, 3] = np.power(10, data[:, 3])
    data[:, 4] = np.power(10, data[:, 4])
    return data

def CalcConfidenceInterval(original, losses, saveDistros = False, percentiles = False):
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
        
def CrossValidateKernel(data):
    bandwidths = np.logspace(-8,-2)
    nfolds = 10
    k_fold = KFold(nfolds)
    score = []
    for bandwidth in bandwidths:
        score.append(0.0)
        for (trainIndices, testIndices) in k_fold.split(data):
            kde = neighbors.KernelDensity(bandwidth = bandwidth)
            kde.fit(data[trainIndices])
            score[-1] = score[-1] + kde.score(data[testIndices])
        score[-1] = -score[-1] / nfolds
        print bandwidth, score[-1]
    plt.plot(bandwidths, score)
    plt.yscale('log')
    plt.show()
    sys.exit()

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

def CleanData(data, minT = 5, maxT = 300, nT = 300, maxT, constG = False, constY = False, constT = False, BoltzmannCorrection = True):
    tlsdata = RemoveOutliers(tlsdata)
    print "Using " + str(tlsdata.shape[0] + " TLS data points."
    if not BoltzmannCorrection:
        tlsdata = RemoveDuplicates(tlsdata)
    nPts = tlsdata.shape[0]
    tMatrix = np.linspace(minT,maxT, num = nT)
    omega = 2 * np.pi * freq
    density = nTLS / vol

    #Import data
    tlsdata[:,5] = 1.0 / tlsdata[:,5] * 0.0101804979 * 1e-12
    tlsdata[:,6] = 1.0 / tlsdata[:,6] * 0.0101804979 * 1e-12
    if constG:
        tlsdata[:,7] = np.mean(tlsdata[:, 7])
        #tlsdata[:, 7] = 1.7
    if constY:
        tlsdata[:, [9,10, 11]] = np.mean(tlsdata[:,[9,10, 11]])
    if constT:
        tlsdata[:,5] = np.mean(tlsdata[:,5])
        tlsdata[:,6] = np.mean(tlsdata[:,6])
    tlsdata[:,[9,10, 11]] = tlsdata[:,[9, 10, 11]] * 6.3227e-7

    return tlsdata
