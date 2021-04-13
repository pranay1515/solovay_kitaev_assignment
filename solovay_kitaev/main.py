import math
import os
import random
from collections import Counter
from random import randint

import numpy as np

from recur_sk import * # will generate the kd-tree as well
from util import dagger

print(
    "\n-------------------SOLOVAY KITAEV ALGORITHM TEST----------------------------\n"
)


def rotateZ(t):
    Z = np.array(
        [
            [math.cos(t / 2) - 1j * math.sin(t / 2), 0],
            [0, math.cos(t / 2) + 1j * math.sin(t / 2)],
        ]
    )
    return Z


def deleteI(str):
    newstr = ""
    for i in range(len(str)):
        if str[i] == "I":
            continue
        newstr = newstr + str[i]
    return newstr


U = rotateZ(randint(1, 100))
V = multiply(solovay_kitaev(U, 0))

print(np.round(U, decimals=3))
print(np.round(V, decimals=3))



numLow = 1 / 3200000000
numHigh = 1 / 3200
listOfNumbers = []

for _ in range(0, 100):
	listOfNumbers.append(random.uniform(numLow, numHigh))

for e in listOfNumbers:
	U = generate_SU2()
	circuit = solovay_kitaev(U, 6)  # check with length 6
	Accuracy = trace_norm(multiply(circuit), U)  # accuracy obtained
	print(
		e,
		Accuracy,
		Counter(deleteI(circuit)),
		sum(Counter(deleteI(circuit)).values()),
	)

