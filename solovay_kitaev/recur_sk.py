# ....................................................MAIN SOLOVAY KITAEV ROUTINE ....................................

import cmath
import math
import pickle
from random import randint

import dill as pickle  # dill is important for dumping lambda functions....
import kdtree
import numpy as np

from approx import multiply
from Gates import *
from Similar import *
from util import dagger

PI = math.acos(-1.0)
basic_gates_sequences = []
l0 = 8  # CRASHES AFTER GATE SEQUENCE OF LENGTH 10


# Creating a class to add a payload to the tuple object to be sent to kdtree
class Sequence(
    object
):  # class to hold the record of the matrix in SU(2) and its payload string
    def __init__(self, r1, c1, r2, c2, payload):
        self.coords = (r1, c1, r2, c2)
        self.payload = payload

    def __len__(self):
        return len(self.coords)

    def __getitem__(self, i):
        return self.coords[i]

    def __repr__(self):
        return "Sequence({},{},{},{},{})".format(
            self.coords[0], self.coords[1], self.coords[2], self.coords[3], self.payload
        )


def GenerateUtil(gates, prefix, n, l):

    if l == 0:
        k = len(basic_gates_sequences)
        # tag the matrix to it
        mat = np.array([[1, 0], [0, 1]])
        print(3 ** 13 - k)  # count down timer
        for i in range(l0):
            if prefix[l0 - i - 1] == "H":
                mat = np.matmul(H, mat)
            if prefix[l0 - i - 1] == "T":
                mat = np.matmul(T, mat)
            if prefix[l0 - i - 1] == "S":
                mat = np.matmul(S, mat)
            # Unroll the matrix and tuple it
        r1 = mat[0][0].real
        c1 = mat[0][0].imag
        r2 = mat[0][1].real
        c2 = mat[0][1].imag

        pts = Sequence(r1, c1, r2, c2, prefix)  # a,bz,by,bx
        basic_gates_sequences.append(pts)

        return

    for i in range(n):
        newPrefix = prefix + gates[i]
        GenerateUtil(gates, newPrefix, n, l - 1)


def Generate_Sequences(gates, l):
    n = len(gates)
    GenerateUtil(gates, "", n, l)


known_gates = [
    "H",
    "I",
    "T",
    "S",
    "s",
    "t",
]  # Hadamard,Identity,Phase-gate and their adjoints...


Generate_Sequences(known_gates, l0)
print(basic_gates_sequences[0])

KD_Tree = kdtree.create(basic_gates_sequences)  # create the KDTree
print(KD_Tree.search_nn([0, 0, 0, 0, 1]))  # nearest neighbour search check!
pin = KD_Tree


def matrix_unroll(Gate):
    r1 = Gate[0][0].real
    c1 = Gate[0][0].imag
    r2 = Gate[0][1].real
    c2 = Gate[0][1].imag
    return [r1, c1, r2, c2]


def basic_approx(U, e0):
    unroll = matrix_unroll(U)
    R = pin.search_nn(unroll)
    coords = R[0].__dict__["data"].__dict__["coords"]
    R_sequence = R[0].__dict__["data"].__dict__["payload"]
    return R_sequence


def dagger_seq(sequence):  # HTHTSST ------> tssthth
    ll = len(sequence)
    new_sequence = ""
    for i in range(ll):
        if sequence[i] == "H":
            new_sequence = "h" + new_sequence
        if sequence[i] == "h":
            new_sequence = "H" + new_sequence
        if sequence[i] == "T":
            new_sequence = "t" + new_sequence
        if sequence[i] == "t":
            new_sequence = "T" + new_sequence
        if sequence[i] == "I":
            new_sequence = "I" + new_sequence
        if sequence[i] == "s":
            new_sequence = "S" + new_sequence
        if sequence[i] == "S":
            new_sequence = "s" + new_sequence
    return new_sequence


e0 = 1 / 32


def solovay_kitaev(U, n):  # input :: gate , depth
    if n == 0:  # base case
        return basic_approx(U, e0)  # returns the sequence
    else:
        U_prev_seq = solovay_kitaev(U, n - 1)  # U(n-1)
        U_prev = multiply(U_prev_seq)
        V, W = decompose(np.matmul(U, dagger(U_prev)))
        V_prev_seq = solovay_kitaev(V, n - 1)
        W_prev_seq = solovay_kitaev(W, n - 1)
        U_now_seq = (
            V_prev_seq
            + W_prev_seq
            + dagger_seq(V_prev_seq)
            + dagger_seq(W_prev_seq)
            + U_prev_seq
        )
        return U_now_seq
