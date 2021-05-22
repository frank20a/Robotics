from sympy import sin, cos, tan, atan, atan2, acos, asin, sqrt, pi, oo


def inverse_kinematic(M06):
    """Given the final 0M6 transition matrix for the working tip returns the angles of the robot joints. The formulas
    are pre-calculated."""
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