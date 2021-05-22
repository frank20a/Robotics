import sympy as sp
from sympy import pi, sqrt

th1, th2, th3, th4, th5, th6 = sp.symbols('th1 th2 th3 th4 th5 th6')
d1 = 183e-3
d4 = 221.5e-3
d6 = 23.7e-3
L2 = 210e-3
L3 = 30e-3
L5 = 5.5e-3
DH_Table = [
    (0,      0,      d1,     th1),
    (-pi/2,  0,      0,      th2-(pi/2)),
    (0,      L2,     0,      th3),
    (-pi/2,  L3,     d4,     th4+pi),
    (pi/2,   0,      0,      th5),
    (-pi/2,  L5,     d6,     th6)
]