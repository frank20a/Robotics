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


def sMf(TransMatrices, s=0, f=6):
    if s < 0 or f > len(TransMatrices): raise IndexError

    res = TransMatrices[s]
    for i in TransMatrices[s+1:f]:
        res = res * i
    return res


def printEval(TransMatrices, theta1, theta2, theta3, theta4, theta5, theta6):
    for i in TransMatrices:
        print(i.evalf(subs={th1: theta1, th2: theta2, th3: theta3, th4: theta4, th5: theta5, th6: theta6}))


def solveDirect(TransMatrices, theta1, theta2, theta3, theta4, theta5, theta6):
    res = np.array([[0], [0], [0], [1]])
    for transMatrix in TransMatrices[::-1]:
        res = transMatrix.evalf(subs={th1: theta1, th2: theta2, th3: theta3, th4: theta4, th5: theta5, th6: theta6}) * res
    return res[0], res[1], res[2]


def solveInverse(M06):
    ((ix, jx, kx, px), (iy, jy, ky, py), (iz, jz, kz, pz), (_, _, _, _)) = M06  # Get data from M06

    # Calculate th1
    if (-237 * ky + 1e4 * py) * (-237 * kx + 1e4 * px) >= 0:
        th1 = atan((-237 * ky + 1e4 * py) / (-237 * kx + 1e4 * px))
    else:
        th1 = -atan(abs((-237 * ky + 1e4 * py) / (-237 * kx + 1e4 * px)))

    # Calculate th3
    A = 21 / 100
    B = px * cos(th1) + py * sin(th1) - (237 / (1e4)) * (kx * cos(th1) + ky * sin(th1))
    C = 183 / 1e3 + (237 * kz / 1e4) - pz
    D = 443 / 2000
    E = 3 / 100
    th3 = asin(((A ** 2 + D ** 2 + E ** 2) - (B ** 2 + C ** 2)) / (2 * A * sqrt(D ** 2 + E ** 2))) + atan(E / D)
    
    # Calculate th2
    F = D - A * sin(th3)
    if C >= 0:
        th2 = asin(F / sqrt(B ** 2 + C ** 2)) - atan(B / C) - th3
    else:
        th2 = -asin(F / sqrt(B ** 2 + C ** 2)) - atan(B / C) - th3

    # Calculate th5
    th5 = -acos(
        kx * cos(th1) * cos(th2) * cos(th3) - kz * cos(th3) * sin(th2) - kz * cos(th2) * sin(th3) + ky * cos(th2) * cos(
            th3) * sin(th1) - kx * cos(th1) * sin(th2) * sin(th3) - ky * sin(th1) * sin(th2) * sin(th3))

    # Calculate th4
    if th5 != 0:
        th4 = asin((ky * cos(th1) - kx * sin(th1)) / sin(th5))
    else:
        G = (py * cos(th1) - px * sin(th1)) - (237 / 1e4) * (ky * cos(th1) - kx * sin(th1))
        if B * sin(th2 + th3) - C * cos(th2 + th3) - A * cos(th3) - E == 0 and G == 0:
            th4 = 0
        else:
            th4 = atan((B * sin(th2 + th3) - C * cos(th2 + th3) - A * cos(th3) - E) / G)

    # Calculate th6
    th6 = asin(
        ix * cos(th4) * sin(th1) - iy * cos(th1) * cos(th4) - iz * cos(th2) * cos(th3) * sin(th4) + iz * sin(th2) * sin(
            th3) * sin(th4) - ix * cos(th1) * cos(th2) * sin(th3) * sin(th4) - ix * cos(th1) * cos(th3) * sin(
            th2) * sin(th4) - iy * cos(th2) * sin(th1) * sin(th3) * sin(th4) - iy * cos(th3) * sin(th1) * sin(
            th2) * sin(th4))

    return th1, th2, th3, th4, th5, th6


def Jacobian(M):
    z0 = sp.Matrix([0, 0, 1])
    z1 = sp.Matrix([M[0][0][2], M[0][1][2], M[0][2][2]])
    t1 = sp.Matrix([M[0][0][3], M[0][1][3], M[0][2][3]])
    z2 = sp.Matrix([sMf(M, 0, 2)[0][2], sMf(M, 0, 2)[1][2], sMf(M, 0, 2)[2][2]])
    t2 = sp.Matrix([sMf(M, 0, 2)[0][3], sMf(M, 0, 2)[1][3], sMf(M, 0, 2)[2][3]])
    z3 = sp.Matrix([sMf(M, 0, 3)[0][2], sMf(M, 0, 3)[1][2], sMf(M, 0, 3)[2][2]])
    t3 = sp.Matrix([sMf(M, 0, 3)[0][3], sMf(M, 0, 3)[1][3], sMf(M, 0, 3)[2][3]])
    z4 = sp.Matrix([sMf(M, 0, 4)[0][2], sMf(M, 0, 4)[1][2], sMf(M, 0, 4)[2][2]])
    t4 = sp.Matrix([sMf(M, 0, 4)[0][3], sMf(M, 0, 4)[1][3], sMf(M, 0, 4)[2][3]])
    z5 = sp.Matrix([sMf(M, 0, 5)[0][2], sMf(M, 0, 5)[1][2], sMf(M, 0, 5)[2][2]])
    t5 = sp.Matrix([sMf(M, 0, 5)[0][3], sMf(M, 0, 5)[1][3], sMf(M, 0, 5)[2][3]])
    # z6 = sp.Matrix([sMf(M, 0, 6)[0][0][2], sMf(M, 0, 6)[0][1][2], sMf(M, 0, 6)[0][2][2]])
    t6 = sp.Matrix([sMf(M, 0, 6)[0][3], sMf(M, 0, 6)[1][3], sMf(M, 0, 6)[2][3]])
    return sp.Matrix([
        [z0.dot(t6)[0],  z1.dot(t6-t1)[0],  z2.dot(t6-t2)[0],  z3.dot(t6-t3)[0],  z4.dot(t6-t4)[0],  z5.dot(t6-t5)[0] ],
        [z0.dot(t6)[1],  z1.dot(t6-t1)[1],  z2.dot(t6-t2)[1],  z3.dot(t6-t3)[1],  z4.dot(t6-t4)[1],  z5.dot(t6-t5)[1] ],
        [z0.dot(t6)[2],  z1.dot(t6-t1)[2],  z2.dot(t6-t2)[2],  z3.dot(t6-t3)[2],  z4.dot(t6-t4)[2],  z5.dot(t6-t5)[2] ],
        [z0[0],          z1[0],             z2[0],             z3[0],             z4[0],             z5[0]            ],
        [z0[1],          z1[1],             z2[1],             z3[1],             z4[1],             z5[1]            ],
        [z0[2],          z1[2],             z2[2],             z3[2],             z4[2],             z5[2]            ]
    ])


M = parseDH2TransMatrix(DH_Table)
print(solveDirect(M, 0, 0, 0, 0, 0, 0))
