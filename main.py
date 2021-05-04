import sympy as sp
import numpy as np
from sympy import sin, cos, tan, atan, atan2, acos, asin, sqrt, pi, oo
from data import *
from pprint import pprint

sp.init_printing(use_unicode=True)


class TransMatrix(sp.Matrix):
    def __init__(self, arg):
        super().__init__()


def DH2TransMatrix(a, L, d, theta):
    return TransMatrix([
        [cos(theta),         -sin(theta),        0,        L        ],
        [sin(theta)*cos(a),  cos(theta)*cos(a),  -sin(a),  -sin(a)*d],
        [sin(theta)*sin(a),  cos(theta)*sin(a),  cos(a),   cos(a)*d ],
        [0,                  0,                  0,        1        ]
    ])


def parseDH2TransMatrix(table):
    res = []
    for a, L, d, theta in table:
        res.append(DH2TransMatrix(a, L, d, theta))
    return tuple(res)


def solveDirect(TransMatrices, theta1, theta2, theta3, theta4, theta5, theta6):
    res = np.array([[0], [0], [0], [1]])
    for transMatrix in TransMatrices[::-1]:
        res = transMatrix.evalf(subs={th1: theta1, th2: theta2, th3: theta3, th4: theta4, th5: theta5, th6: theta6}) * res
    return res[0], res[1], res[2]


def printEval(TransMatrices, theta1, theta2, theta3, theta4, theta5, theta6):
    for i in TransMatrices:
        print(i.evalf(subs={th1: theta1, th2: theta2, th3: theta3, th4: theta4, th5: theta5, th6: theta6}))


def sMf(TransMatrices, s=0, f=6):
    if s < 0 or f > len(TransMatrices): raise IndexError

    res = TransMatrices[s]
    for i in TransMatrices[s+1:f]:
        res = res * i
    return res


def generateEquations(M):
    res = []
    for i in range(0, len(M)-1):
        for j in range(i, -1, -1):
            pass


M = parseDH2TransMatrix(DH_Table)
print(solveDirect(M, 0, 0, 0, 0, 0, 0))

pprint(M[0]**-1 * sMf(M, 0, 6)-sMf(M, 1, 6))
pprint(M[4]**-1 * M[3]**-1 * M[2]**-1 * M[1]**-1 * M[0]**-1 * sMf(M, 0, 6)-sMf(M, 5, 6))
