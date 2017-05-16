import numpy as np
import sys

def ImportData(fn):
    inStream = open(fn, 'r')
    data = []
    while True:
        line = inStream.readline()
        if not line:
            break
        row = line.split()
        if len(row) == 12:
            data.append(map(float, row))
    return np.array(data)


def RemoveDuplicates(mat, cut):
    colsToCompare = [0, 1, 2]
    cutMat = np.copy(mat)
    nCut = 0
    nRows = mat.shape[0]
    i = 0
    while i < cutMat.shape[0]:
        j = 0
        while j < i:
            sqDiff = 0.0
            for ci in colsToCompare:
                if cutMat[j - nCut, ci] < 1e8:
                    sqDiff = sqDiff + (cutMat[j - nCut, ci] - cutMat[i, ci]) ** 2
                else:
                    sqDiff = sqDiff + (cutMat[j - nCut, ci] - cutMat[i, ci]) ** 2 / cutMat[j - nCut, ci] ** 2
            diff = np.sqrt(sqDiff)
            if(diff < cut):
                cutMat = np.delete(cutMat, i, 0)
                nCut = nCut + 1
                break
            j = j + 1
        i = i + 1
    print nCut
    return cutMat

if len(sys.argv) > 1:
    suff = "_" + sys.argv[1]
else:
    suff = ""

#tlsData = np.loadtxt("tls.tot")
tlsData = ImportData("tls.tot" + suff)
#tlsData = RemoveDuplicates(tlsData, 0.001)
np.savetxt("pp-tls.tot" + suff, tlsData)
