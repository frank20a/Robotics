import sympy as sp
from sympy import pi, sqrt

th1, th2, th3, th4, th5, th6 = sp.symbols('th1 th2 th3 th4 th5 th6')
d1 = 0
d4 = 0
d5 = 0
L1 = 0
L2 = 0
L5 = 0
DH_Table = [
    (0,     0,     d1,    th1),
    (-pi/4, L1,    0,     th2),
    (0,     L2,    0,     th3),
    (-pi/4, 0,     d4,    th4),
    (0,     0,     d5,    th5),
    (pi/4,  L5,    0,     th6)
]