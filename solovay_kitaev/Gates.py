#   ----------------------------------------  Quantum Gates ---------------------------------------------------------------------

import math

import numpy as np

from phase_balance import *

PI = math.acos(-1.0)

# Library Gates ............................
H = np.array([[1, 1], [1, -1]]) / math.sqrt(2)  # Hadamard
T = np.array([[1, 0], [0, math.cos(PI / 4) + 1j * math.sin(PI / 4)]])  # Phase Gate
I = np.array([[1, 0], [0, 1]])  # Identity matrix
S = np.matmul(T, T)
s = S.transpose().conjugate()  # S dagger as s
t = T.transpose().conjugate()  # T dagger as t
h = H.transpose().conjugate()


# PAULI GATES
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])


# Convert all the unitary matrix to a matrix in SU(2)
H = SU2(H)
h = SU2(h)
T = SU2(T)
I = SU2(I)
S = SU2(S)
s = SU2(s)
t = SU2(t)
X = SU2(X)
Y = SU2(Y)
Z = SU2(Z)
