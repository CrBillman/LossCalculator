import numpy as np
import matplotlib.pyplot as plt

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
        return 0
    else:
        coshTerm = np.cosh(np.sqrt(delta**2 + delta0**2)/(2 * kb * T))
    return 1/coshTerm

def dipoleMoment(delta, delta0, gamma):
    return delta**2 / (delta**2 + delta0**2) * gamma **2

def prefactor(Modulus, T):
    kb = 8.617342E-5
    return 1/(3 * kb * T * Modulus)

def calculateLossFunction(data, density, omega, tMatrix):
    qMatrix = np.zeros(tMatrix.shape)
    nTemp = tMatrix.shape[0]
    try:
        nPts = data.shape[0]
    except AttributeError:
        nPts = len(data)
    
    norm = density / float(nPts)
    for it in range(0, nTemp):
        for ip in range(0,len(data)):
            delta = data[ip, 0]
            V = data[ip, 1]
            tau1 = data[ip, 5]
            tau2 = data[ip, 6]
            Y = data[ip, 10]
            gamma = data[ip, 7]
            qMatrix[it] = qMatrix[it] + prefactor(Y,tMatrix[it]) * resonanceFraction(tau1, tau2, V, tMatrix[it], delta, omega) * sechTerm(delta, 1e-4, tMatrix[it]) * dipoleMoment(delta, 1e-4, gamma)
    qMatrix[:] = qMatrix[:] * norm * 1e3
    return qMatrix

def RemoveOutliers(data):
    #q = data[np.percentile(data[:, 7] , 99), 7]
    qL = np.percentile(data[:, 7], 2)
    qU = np.percentile(data[:, 7], 98)
    i = 0
    print(data.shape)
    while i < data.shape[0]:
        if data[i, 7] < qL:
            data = np.delete(data, (i), axis = 0)
        if data[i, 7] > qU:
            data = np.delete(data, (i), axis = 0)
        else:
            i = i + 1
    print(data.shape)
    return data

def RemoveDuplicates(data):
    window = 200
    colsToComp = [0, 2, 3, 4]
    eComp = 0.001
    dComp = 0.05
    #colsToComp = [0, 1]
    cut = 0.00
    absCut = {0 : eComp, 
            2 : dComp, 
            3 : eComp,
            4 : eComp}
    nRows = data.shape[0]
    noDups = []
    for i in range(0, nRows):
        match = False
        #for j in range(i + 1, nRows):
        ul = min(i + window, nRows)
        for j in range(i + 1, ul):
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
