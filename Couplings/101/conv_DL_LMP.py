import sys
import numpy as np
#import scipy.linalg

def ConvertDL2LMP(fnDL, fnLMP, step):
    dl = open(fnDL, 'r')
    lmp = open(fnLMP, 'wb')

    lmp.write("ITEM: TIMESTEP\n" + step + '\n')

    dl.readline()
    wRows, cType, nAtoms = dl.readline().split()
    nAtoms = "1008"
    lmp.write("ITEM: NUMBER OF ATOMS\n" + nAtoms + '\n')

    cVecs = np.empty([3,3])
    for i in xrange(0, 3):
        row = dl.readline().split()
        for j in xrange(0,3):
            cVecs[i, j] = float(row[j])
    print cVecs
    #hMat = scipy.linalg.lu(cVecs)[2]
    #hMat, rot = scipy.linalg.schur(cVecs, output = 'complex')
    #print hMat, rot
    rot, hMat = LUDecomp(cVecs)
    lmp.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
    lmp.write("0.0 " + str(hMat[0,0]) + " " + str(hMat[0,1]) + '\n')
    lmp.write("0.0 " + str(hMat[1,1]) + " " + str(hMat[0,2]) + '\n')
    lmp.write("0.0 " + str(hMat[2,2]) + " " + str(hMat[1,2]) + '\n')

    lmp.write("ITEM: ATOMS id type x y z\n")
    aDict = {}
    aTypes = 0
    while True:
        line = dl.readline()
        if not line:
            break
        aName, aNum = line.split()
        try:
            aIndex = aDict[aName]
        except KeyError:
            aDict[aName] = str(aTypes+1)
            aIndex = aDict[aName]
            aTypes = aTypes + 1
        x, y, z = dl.readline().split()
        RotateCoordinates(x, y, z, hMat, rot)
        for i in xrange(0, int(wRows)):
            dl.readline()
        lmp.write(" ".join([aNum, aIndex, x, y, z, '\n']))

    lmp.close()
    dl.close()
    print aTypes
    return str(aTypes)

def ConvertDL2Data(fnDL, fnLMP, aTypes):
    dl = open(fnDL, 'r')
    lmp = open(fnLMP, 'wb')

    lmp.write(dl.readline())
    lmp.write('\n')

    wRows, cType, nAtoms = dl.readline().split()
    nAtoms = "1008"
    lmp.write(nAtoms + " atoms\n")
    lmp.write(aTypes + " atom types\n")
    lmp.write('\n')

    cVecs = np.empty([3,3])
    for i in xrange(0, 3):
        row = dl.readline().split()
        for j in xrange(0,3):
            cVecs[i, j] = float(row[j])
    #hMat = scipy.linalg.lu(cVecs)[2]
    rot, hMat = LUDecomp(cVecs)
    lmp.write("0.0 " + str(hMat[0,0]) + " xlo xhi\n")
    lmp.write("0.0 " + str(hMat[1,1]) + " ylo yhi\n")
    lmp.write("0.0 " + str(hMat[2,2]) + " zlo zhi\n")
    lmp.write(" ".join([str(hMat[0,1]), str(hMat[0,2]), str(hMat[1,2]), "xy", "xz", "yz", '\n']))
    lmp.write('\n')

    charges, masses = ReadField()

    lmp.write("Masses\n\n")
    for i in xrange(0, len(masses)):
        lmp.write(str(i+1) + ' ' + masses[i] + '\n')
    lmp.write('\n')

    lmp.write("Atoms # charge\n\n")
    aDict = {}
    aTypes = 0
    while True:
        line = dl.readline()
        if not line:
            break
        aName, aNum = line.split()
        try:
            aIndex = aDict[aName]
            aCharge = charges[int(aIndex) - 1]
        except KeyError:
            aDict[aName] = str(aTypes+1)
            aIndex = aDict[aName]
            aCharge = charges[int(aIndex) - 1]
            aTypes = aTypes + 1
        x, y, z = dl.readline().split()
        RotateCoordinates(x, y, z, rot, hMat)
        for i in xrange(0, int(wRows)):
            dl.readline()
        lmp.write(" ".join([aNum, aIndex, aCharge, x, y, z, '\n']))


    return

def ReadField():
    try:
        fd = open("FIELD", 'r')
        charges = []
        masses = []
        while True:
            line = fd.readline()
            if not line:
                break
            if "atoms" in line:
                sLine = fd.readline().split()
                masses.append(sLine[1])
                charges.append(sLine[2])
    except IOError:
        print "There is no FIELD file, using the masses and charged defined in the python code."
        charges = ["-1.2", "2.5"]
        masses = ["180.95000", "15.99940"]
    return charges, masses

def LUDecomp(matrix):
    u = np.zeros(matrix.shape)
    l = np.zeros(matrix.shape)
    for i in xrange(0, matrix.shape[0]):
        #for j in xrange(0, matrix.shape[1]):
        for j in xrange(i, matrix.shape[0]):
            uSum = 0
            for k in xrange(0, i -1):
                uSum = uSum + u[k, j] * l[i, k]
            u[i, j] = matrix[i, j] - uSum
        for j in xrange(0, i+1):
            lSum = 0
            for k in xrange(0, i -1):
                lSum = lSum + u[k, j] * l[i, k]
            l[i, j] = (matrix[i, j] - lSum) / u[j, j]
    return l, u

def RotateCoordinates(x, y, z, C, R):
    floatList = map(float, [x, y, z])
    xVec = np.transpose(np.matrix(floatList))
    #xVec = np.matrix(floatList)
    #xpVec = np.matmul(R, xVec)
    r1 = np.matmul(np.linalg.inv(R), C)
    r2 = np.matmul(R, np.linalg.inv(C))
    rot = np.matmul(r1, r2)
    xpVec = np.matmul(rot, xVec)
    strList = map(str, list(xpVec))
    return xpVec
    #xVec = np.empty([3])
    #xVec[0] = x
    #xVec[1] = y
    #xVec[2] = float(z)


try:
    ifn = sys.argv[1]
except IndexError:
    print "Need to specify a DL_POLY CONFIG file."
    sys.exit(-1)

try:
    ofn = sys.argv[2]
except IndexError:
    print "Need to specify a LAMMPS file to write to."
    sys.exit(-1)
try:
    step = sys.argv[3]
except IndexError:
    step = '0'


aTypes = ConvertDL2LMP(ifn, ofn + ".dump", step)
ConvertDL2Data(ifn, ofn + ".dat", aTypes)
