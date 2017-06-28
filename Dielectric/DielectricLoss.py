import numpy as np

def Tau(Tau1, Tau2, V, T, delta):
    numCut = 80
    kb = 8.617342E-5
    if (V+delta) / (kb * T)  > numCut or (V - delta) / (kb * T) > numCut:
        return np.exp(numCut)
    else:
        return 0.5 * (Tau1 * np.exp((V+delta)/(kb * T)) + Tau2 * np.exp((V-delta)/(kb * T)))

def suscResonance(Tau1, Tau2, V, T, delta, omega):
    tau = Tau(Tau1, Tau2, V, T, delta)
    if(np.isnan(tau)):
        return 0.0
    else:
        return 1.0 / (1.0 - omega * tau * 1j)    

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

def suscPrefactor(T):
    kb = 8.617342E-5
    denom = 3.0 * kb * T
    denom = denom * 1.60218e-19
    return 1.0 / (denom)

def calculateSusceptibility(data, density, temp, fMatrix, BoltzmannCorrection = False):
    kb = 8.617342E-5
    qMatrix = np.zeros(fMatrix.shape, dtype=complex)
    nFreq = fMatrix.shape[0]
    try:
        nPts = data.shape[0]
    except AttributeError:
        nPts = len(data)
    if BoltzmannCorrection:
        BCc = 0
        norm = 0.0
        for ip in xrange(0, len(data)):
            V = data[ip, 1]
            norm = norm + np.exp( V / (kb * 1000))
        #norm = density / (norm)
        norm = density / float(nPts)
        print norm, density / float(nPts)
    else:
        BCc = 1
        norm = density / float(nPts)
    print norm

    for fi in xrange(0, nFreq):
        for ip in xrange(0,len(data)):
            delta = data[ip, 0]
            V = data[ip, 1]
            tau1 = data[ip, 5]
            tau2 = data[ip, 6]
            Y = data[ip, 11]
            gamma = data[ip, 7]
            omega = fMatrix[fi]
            BC = max(BCc, np.exp(V / (kb * 1000)) )
            qMatrix[fi] = qMatrix[fi] + BC * suscPrefactor(temp) * suscResonance(tau1, tau2, V, temp, delta, omega) * sechTerm(delta, 1e-4, temp) * dipoleMoment(delta, 1e-4, gamma)
    qMatrix[:] = qMatrix[:] * norm
    return qMatrix
