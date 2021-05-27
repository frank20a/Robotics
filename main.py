import numpy as np
from sympy import sin, cos, tan, atan, atan2, acos, asin, sqrt, pi, oo, sign
from data import *
import pandas as pd

sp.init_printing(use_unicode=True)


class TransMatrix(sp.Matrix):
    """TransMatrix subclasses sympy's Matrix in order to add extra features. Currently not implemented!"""

    def __init__(self, arg):
        super().__init__()


def DH2TransMatrix(a, L, d, theta):
    """Creates Transition Matrix for a specific joint given its DH-Parameters"""
    return TransMatrix([
        [cos(theta), -sin(theta), 0, L],
        [sin(theta) * cos(a), cos(theta) * cos(a), -sin(a), -sin(a) * d],
        [sin(theta) * sin(a), cos(theta) * sin(a), cos(a), cos(a) * d],
        [0, 0, 0, 1]
    ])


def parseDH2TransMatrix(table):
    """Creates list of Transition Matrices given the DH-Table for the robot. list[0] is 0M1, list 1 is 1M2 etc."""
    res = []
    for a, L, d, theta in table:
        res.append(DH2TransMatrix(a, L, d, theta))  # Call DH2TransMatrix recursively
    return tuple(res)


def sMf(TransMatrices, s: int = 0, f: int = 6):
    """Returns the Transition Matrix from joint s to joint f (sMs+1 * s+1Ms+2 * ... * f-1Mf)."""
    if s < 0 or f > len(TransMatrices): raise IndexError

    res = TransMatrices[s]
    for i in TransMatrices[s + 1:f]:
        res = res * i   # Multiply TransMatrix's
    return res


def printEval(TransMatrices, theta: tuple):
    """Print the Transition Matrices after evaluating symbolic variables with constants."""
    theta1, theta2, theta3, theta4, theta5, theta6 = theta

    for i in TransMatrices:
        print(i.evalf(subs={th1: theta1, th2: theta2, th3: theta3, th4: theta4, th5: theta5, th6: theta6}))


def solveDirect(TransMatrices, theta):
    """Given the list of Transition Matrices and the joint angles, return the position of the working tip."""
    theta1, theta2, theta3, theta4, theta5, theta6 = theta

    res = np.array([
        [1, 1, 1, 0],
        [1, 1, 1, 0],
        [1, 1, 1, 0],
        [0, 0, 0, 1]
    ])

    for transMatrix in TransMatrices[::-1]:
        res = np.array(transMatrix.evalf(
            subs={th1: theta1, th2: theta2, th3: theta3, th4: theta4, th5: theta5, th6: theta6})) * res

    return res


def inverse_kinematic(M06, TransMatrices):
    """Given the final 0M6 transition matrix for the working tip returns the angles of the robot joints. The formulas
    are pre-calculated."""

    solutions = []
    tests = inverse_solutions(M06, TransMatrices)

    for test in tests:
        test = [complex(v) for v in test]

        if any(map(lambda x: x.imag != 0, test)):
            continue
        else:
            test = [v.real for v in test]

        error = sum(sum(abs(M06 - solveDirect(TransMatrices, test)))) > 1e-2
        lims = bool(any((abs(test[0]) > 175*pi/180,
                   -pi/2 > test[1] > 36.7*pi/180,
                    -80*pi/180 > test[2] > pi/2,
                   test[3] > 175*pi/180,
                   test[4] < -100*pi/180,
                   abs(test[5]) > 147.5*pi/180)))

        if not lims:
            solutions.append(test)

    # print(solutions)
    return np.array(solutions[0])


def inverse_solutions(M06, TransMatrices):
    """Solves inverse kinematic problem and returns list with all the possible solutions"""
    ((ix, jx, kx, px), (iy, jy, ky, py), (iz, jz, kz, pz), (_, _, _, _)) = M06  # Get data from M06

    # Calculate 2 th1 solutions
    th1s = []
    tmp = atan((-237 * ky + 1e4 * py) / (-237 * kx + 1e4 * px))
    th1s.append([tmp])
    if (-237 * ky + 1e4 * py) / (-237 * kx + 1e4 * px) > 0:
        th1s.append([tmp - pi])
    elif (-237 * ky + 1e4 * py) / (-237 * kx + 1e4 * px) < 0:
        th1s.append([tmp + pi])

    # Calculate 2 th3 solutions (1 per)
    A = 21 / 100
    C = 183 / 1e3 + (237 * kz / 1e4) - pz
    D = 443 / 2000
    E = 3 / 100

    th3s = []
    for th1, in th1s:
        B = px * cos(th1) + py * sin(th1) - (237 / 1e4) * (kx * cos(th1) + ky * sin(th1))
        th3s.append([th1,
                     asin(((A ** 2 + D ** 2 + E ** 2) - (B ** 2 + C ** 2)) / (2 * A * sqrt(D ** 2 + E ** 2))) + atan(
                         E / D)])

    # Calculate 4 th2 solutions (2 per)
    th2s = []
    for th1, th3 in th3s:
        B = px * cos(th1) + py * sin(th1) - (237 / 1e4) * (kx * cos(th1) + ky * sin(th1))
        F = D - A * sin(th3)

        if C != 0 and B != 0:
            tmp1 = sign(C) * asin(F / sqrt(B ** 2 + C ** 2)) - sign(B) * sign(C) * atan(abs(B / C)) - th3
            tmp2 = (pi - sign(C) * asin(F / sqrt(B ** 2 + C ** 2))) - sign(B) * sign(C) * atan(abs(B / C)) - th3
        elif C == 0:
            tmp1 = acos(F / B) - th3
            tmp2 = -acos(F / B) - th3
        elif B == 0:
            tmp1 = asin(F / C) - th3
            tmp2 = (pi - asin(F / C)) - th3
        th2s.append([th1, tmp1, th3])
        th2s.append([th1, tmp2, th3])

    # Calculate 8 th5 solutions (2 per)
    th5s = []
    for th1, th2, th3 in th2s:
        tmp = acos(kx * cos(th1) * cos(th2) * cos(th3) - kz * cos(th3) * sin(th2) - kz * cos(th2) * sin(th3) + ky * cos(
            th2) * cos(th3) * sin(th1) - kx * cos(th1) * sin(th2) * sin(th3) - ky * sin(th1) * sin(th2) * sin(th3))
        th5s.append([th1, th2, th3, tmp])
        th5s.append([th1, th2, th3, -tmp])

    # Calculate 16 th4 solutions (2 per)
    th4s = []
    for th1, th2, th3, th5 in th5s:
        if th5 != 0:
            tmp = asin((ky * cos(th1) - kx * sin(th1)) / sin(th5))
            th4s.append([th1, th2, th3, tmp, th5])
            th4s.append([th1, th2, th3, pi - tmp, th5])
        else:
            th4s.append([th1, th2, th3, 0, th5])

    # Calculate 32 th6 solutions (2 per)
    solutions = []
    for th1, th2, th3, th4, th5 in th4s:
        tmp = asin(ix * cos(th4) * sin(th1) - iy * cos(th1) * cos(th4) - iz * cos(th2) * cos(
            th3) * sin(th4) + iz * sin(th2) * sin(th3) * sin(th4) - ix * cos(th1) * cos(
            th2) * sin(th3) * sin(th4) - ix * cos(th1) * cos(th3) * sin(th2) * sin(
            th4) - iy * cos(th2) * sin(th1) * sin(th3) * sin(th4) - iy * cos(th3) * sin(
            th1) * sin(th2) * sin(th4))

        solutions.append((th1, th2, th3, th4, th5, tmp))
        solutions.append((th1, th2, th3, th4, th5, pi - tmp))

    return solutions


def Jacobian(M):
    """Given the list of the Transition Matrices returns the Jacobian Matrix for the robot."""

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
        [z0.dot(t6)[0], z1.dot(t6 - t1)[0], z2.dot(t6 - t2)[0], z3.dot(t6 - t3)[0], z4.dot(t6 - t4)[0],
         z5.dot(t6 - t5)[0]],
        [z0.dot(t6)[1], z1.dot(t6 - t1)[1], z2.dot(t6 - t2)[1], z3.dot(t6 - t3)[1], z4.dot(t6 - t4)[1],
         z5.dot(t6 - t5)[1]],
        [z0.dot(t6)[2], z1.dot(t6 - t1)[2], z2.dot(t6 - t2)[2], z3.dot(t6 - t3)[2], z4.dot(t6 - t4)[2],
         z5.dot(t6 - t5)[2]],
        [z0[0], z1[0], z2[0], z3[0], z4[0], z5[0]],
        [z0[1], z1[1], z2[1], z3[1], z4[1], z5[1]],
        [z0[2], z1[2], z2[2], z3[2], z4[2], z5[2]]
    ])


def getAnglePolynomial(th0, thf, tf):
    """Returns a lambda that implements a 3rd degree polynomial for the angle trajectory provided the trajectory
    starting and ending angle and the time duration"""
    return lambda t: th0 + (3 * (thf - th0) / (tf**2)) * (t**2) + (2 * (th0 - thf) / (tf**3)) * (t**3)


def get2PTrajectory(TransMatrices, start, finish, T, Ts=0.05, filename=None):
    t = np.arange(0, T+Ts, Ts)

    s = inverse_kinematic(start, TransMatrices)
    print(s)
    f = inverse_kinematic(finish, TransMatrices)

    polynomials = []
    for th0, thf in zip(s, f): polynomials.append(getAnglePolynomial(th0, thf, T))

    trajectory = pd.DataFrame(columns=['th1', 'th2', 'th3', 'th4', 'th5', 'th6'], dtype='float64')
    # trajectory.t = t
    for col, polynomial in zip(trajectory.columns, polynomials):
        trajectory[col] = polynomial(t)

    if filename is not None:
        with open(str('./trajectories/' + filename + '.csv'), 'w', encoding='utf-8') as file:
            trajectory.to_csv(path_or_buf=file, index=False)
    else:
        return trajectory


def getTrajectoryCoeffs(s, p1, p2, f, t):
    t2, t3, t4 = t

    # Calculate angular speeds at intermediate points
    dP1 = ((6 * (t3 + t4) / (t2 * t3)) * ((t2 ** 2) * (p2 - p1) + (t3 ** 2) * (p1 - s)) - (3 * t2 / (t3 * t4)) * (
            (t3 ** 2) * (f - p2) + (t4 ** 2) * (p2 - p1))) / (
                  4 * t2 * t3 + 3 * t2 * t4 + 4 * t3 * t4 + 4 * t3 ** 2)

    dP2 = (-(3 * t4 / (t2 * t3)) * ((t2 ** 2) * (p2 - p1) + (t3 ** 2) * (p1 - s)) + (6 * (t2 + t3) / (t3 * t4)) * (
            (t3 ** 2) * (f - p2) + (t4 ** 2) * (p2 - p1))) / (
                  4 * t2 * t3 + 3 * t2 * t4 + 4 * t3 * t4 + 4 * t3 ** 2)

    dF = np.zeros(6)

    # Calculate Polynomial Coeffs
    a = np.zeros((3, 6))
    b = np.zeros((3, 6))
    c = np.zeros((3, 6))
    d = np.zeros((3, 6))

    a[0, :] = s
    c[0, :] = (3 / (t2 ** 2)) * (p1 - s) - (1 / t2) * dP2
    d[0, :] = (-2 / (t2 ** 3)) * (p1 - s) + (1 / (t2 ** 2)) * dP2

    a[1, :] = p1
    b[1, :] = dP2
    c[1, :] = (3 / (t3 ** 2)) * (p2 - p1) - (2 / t3) * dP1 - (1 / t3) * dP2
    d[1, :] = (-2 / (t3 ** 3)) * (p2 - p1) + (1 / (t3 ** 2)) * (dP2 + dP1)

    a[2, :] = p2
    b[2, :] = dP2
    c[2, :] = (3 / (t4 ** 2)) * (f - p2) - (2 / t4) * dP2 - (1 / t4) * dF
    d[2, :] = (-2 / (t4 ** 3)) * (f - p2) + (1 / (t4 ** 2)) * (dF + dP2)

    return a, b, c, d


def getTrajectory(TransMatrices, S, P1, P2, F, t: tuple, Ts=0.05, degrees: bool = False, filename: str = None,
                  humane: bool = False):

    t2, t3, t4 = t
    T12 = np.arange(0, t2, Ts)
    T23 = np.arange(t2, t3, Ts)
    T34 = np.arange(t3, t4, Ts)

    # Solve Inverse Kinematic for trajectory points
    s = inverse_kinematic(S, TransMatrices)
    p1 = inverse_kinematic(P1, TransMatrices)
    p2 = inverse_kinematic(P2, TransMatrices)
    f = inverse_kinematic(F, TransMatrices)

    # Get Coeffs
    A, B, C, D = getTrajectoryCoeffs(s, p1, p2, f, t)

    # Create polynomials
    trajectory = pd.DataFrame(columns=['th1', 'th2', 'th3', 'th4', 'th5', 'th6'], dtype='float64')

    for i, Tt in enumerate((T12, T23, T34)):
        temp = pd.DataFrame(columns=['th1', 'th2', 'th3', 'th4', 'th5', 'th6'], dtype='float64')
        for n, col in enumerate(temp.columns):
            temp[col] = (A[i, n] + B[i, n]*Tt + C[i, n]*(Tt**2) + D[i, n]*(Tt**3)) * (180 if degrees else 1)
        trajectory = trajectory.append(temp, ignore_index=True)

    # print(trajectory)

    if filename is not None:
        with open(str('./trajectories/' + filename + '.csv'), 'w', encoding='utf-8') as file:
            trajectory.to_csv(path_or_buf=file, index=False, header=humane, line_terminator='\n')
    else:
        return trajectory


if __name__ == '__main__':
    M = parseDH2TransMatrix(DH_Table)

    Start = np.array([
        [0,  0,  1,  (287.5 - 83.62) * 1e-3],
        [0, -1,  0,  62.5e-3],
        [1,  0,  0,  275e-3],
        [0,  0,  0,  1]
    ])

    P1 = np.array([
        [0,  0,  1,  (287.5 - 83.62) * 1e-3],
        [0, -1,  0,  62.5e-3],
        [1,  0,  0,  300e-3],
        [0,  0,  0,  1]
    ])

    P2 = np.array([
        [0,  0,  1,  (170 - 83.62) * 1e-3],
        [-1, 0,  0,  62.5e-3],
        [0, -1,  0,  300e-3],
        [0,  0,  0,  1]
    ])

    Finish = np.array([
        [0,  0, 1,  (170 - 83.62) * 1e-3],
        [-1, 0, 0,  220e-3],
        [0, -1, 0,  300e-3],
        [0,  0, 0,  1]
    ])

    getTrajectory(M, Start, P1, P2, Finish, (2, 4, 6), filename='test1')
