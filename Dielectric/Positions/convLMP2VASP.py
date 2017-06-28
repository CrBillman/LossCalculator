import numpy as np
import sys

ifn = sys.argv[1]
ofn = sys.argv[2]

atomTypes = ['Si', 'O']

lmp = open(ifn, 'r')

lmp.readline()
timestep = int(lmp.readline())
lmp.readline()
nAtoms = int(lmp.readline())
lmp.readline()
hMatrix = np.zeros([3,3])
fsLine = map(float, lmp.readline().split())
hMatrix[0, 0] = fsLine[1] - fsLine[0]
hMatrix[0, 1] = fsLine[2]
fsLine = map(float, lmp.readline().split())
hMatrix[1, 1] = fsLine[1] - fsLine[0]
hMatrix[0, 2] = fsLine[2]
fsLine = map(float, lmp.readline().split())
hMatrix[2, 2] = fsLine[1] - fsLine[0]
hMatrix[1, 2] = fsLine[2]

lmp.readline()
pos = np.empty([nAtoms, 4])
atomNumbers = []
for i in xrange(0, nAtoms):
    sfLine = map(float, lmp.readline().split())
    pos[i, 1:] = sfLine[2:]
    atomIndex = int(sfLine[1]) - 1
    pos[i, 0] = atomIndex
    try:
        atomNumbers[atomIndex] = atomNumbers[atomIndex] + 1
    except IndexError:
        atomNumbers.append(1)

pos = pos[np.argsort(pos[:, 0])]

lmp.close()

vasp = open(ofn, 'w')
vasp.write("Converted from " + ifn + '\n')
vasp.write("1.0\n")
for dim in xrange(0, 3):
    vasp.write("  ".join(map('{:10.8f}'.format, hMatrix[dim, :])) + '\n')

vasp.write(" ".join(atomTypes) + '\n')
vasp.write(" ".join(map(str, atomNumbers)) + '\n')
vasp.write("Cartesian\n")
for i in xrange(0, nAtoms):
    vasp.write(" ".join(map('{:10.8f}'.format, pos[i, 1:])) + '\n')
vasp.close()
