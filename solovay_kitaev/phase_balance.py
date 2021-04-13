# 	GLOBAL PHASE BALANCER
# 	U  be any unitary matrix st UU*=I .
# 	U' be its transform in SU(2).
# 	then U' = [1/det(U)]^0.5

import cmath
import math

import numpy as np

from util import determinant


def SU2(U):
    t = complex(0, 0)
    t = t + determinant(U)
    globalPhase = (1 / t) ** 0.5
    return U * globalPhase
