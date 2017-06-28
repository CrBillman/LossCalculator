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
        BC = np.empty([nPts])
        norm = 0.0
        for ip in xrange(0, nPts):
            V1 = data[ip, 3]
            norm = norm + np.exp( (V1) / (kb * 1000))
            BC[ip] = np.exp( (V1) / (kb * 1000))
        #norm = density / (norm)
        norm = density / float(nPts)
    else:
        BC = np.ones([nPts])
        norm = density / float(nPts)
    
    for it in xrange(0, nTemp):
        for ip in xrange(0,nPts):
            delta = data[ip, 0]
            #V = data[ip, 1]
            V = 0.5 * (data[ip, 3] + data[ip, 4])
            V1 = data[ip, 3]
            tau1 = data[ip, 5]
            tau2 = data[ip, 6]
            Y = data[ip, 11]
            gamma = data[ip, 7]
            qMatrix[it] = qMatrix[it] + BC[ip] * prefactor(Y,tMatrix[it]) * resonanceFraction(tau1, tau2, V, tMatrix[it], delta, omega) * sechTerm(delta, 1e-4, tMatrix[it]) * dipoleMoment(delta, 1e-4, gamma)
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

def plotLoss(tMatrix, originalLoss, lLoss, uLoss):
    plt.plot(tMatrix, originalLoss, color = 'blue', label = "Calculated Loss")
    plt.plot(tMatrix, lLoss, color = 'orange', label = "95% Confidence Interval", lw = 2, ls = '--')
    plt.plot(tMatrix, uLoss, color = 'orange', lw = 2, ls = '--')
    plt.legend()
    plt.ylim([0.0, np.amax(uLoss)])
    plt.show()

def saveLoss(tMatrix, originalLoss, lLoss, uLoss, fn):
    prntMat = np.empty([tMatrix.shape[0], 4])
    prntMat[:, 0] = tMatrix[:]
    prntMat[:, 1] = originalLoss[:]
    prntMat[:, 2] = lLoss[:]
    prntMat[:, 3] = uLoss[:]
    np.savetxt(fn, prntMat)
