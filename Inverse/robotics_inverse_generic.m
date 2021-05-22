clear;
clc;

% ORIENTATION AND DISPLACEMENT PARAMETERS 

% ORIENTATION 
ix = 0;
iy = 0;
iz = 1;

jx = 0;
jy = -1;
jz = 0;

kx = 1;
ky = 0;
kz = 0;

% position

px = 0.3;
py = 0.3;
pz = 0.3;

% POSITION AND ORIENTATION OF THE WORK EDGE (SYSTEM X6Y6Z6)

M06 = [ix jx kx px ; iy jy ky py ; iz jz kz pz ; 0 0 0 1];

%% THETA1 - 2 SOLUTIONS

% CALCULATION OF THETA1 (SOLUTION 1)

th1_sol1 = atan((-237*ky+(1e4)*py)/(-237*kx+(1e4)*px));

% CALCULATION OF THETA1 (SOLUTION 2)

if ((-237*ky+(1e4)*py)/(-237*kx+(1e4)*px) > 0)
    th1_sol2 = atan((-237*ky+(1e4)*py)/(-237*kx+(1e4)*px)) - pi;
elseif ((-237*ky+(1e4)*py)/(-237*kx+(1e4)*px) < 0)
    th1_sol2 = atan((-237*ky+(1e4)*py)/(-237*kx+(1e4)*px)) + pi;
else
    th1_sol2 = th1_sol1;
end

%% THETA3 - 2 SOLUTIONS

% CALCULATION OF THETA3 (SOLUTION 1)

A = 21/100;
B = px*cos(th1_sol1) + py*sin(th1_sol1) - (237/(1e4))*(kx*cos(th1_sol1) + ky*sin(th1_sol1));
C = 183/(1e3) + (237*kz/(1e4)) - pz;

D = 443/2000;
E = 3/100;

th3_sol1 = asin(((A^2 + D^2 + E^2) - (B^2 + C^2))/(2*A*sqrt(D^2 + E^2)))+atan(E/D);

% CALCULATION OF THETA3 (SOLUTION 2)

A = 21/100;
B = px*cos(th1_sol2) + py*sin(th1_sol2) - (237/(1e4))*(kx*cos(th1_sol2) + ky*sin(th1_sol2));
C = 183/(1e3) + (237*kz/(1e4)) - pz;

D = 443/2000;
E = 3/100;

th3_sol2 = asin(((A^2 + D^2 + E^2) - (B^2 + C^2))/(2*A*sqrt(D^2 + E^2)))+atan(E/D);

%% THETA2 - 4 SOLUTIONS 

% CALCULATION OF THETA2 (SOLUTION 1)

B = px*cos(th1_sol1) + py*sin(th1_sol1) - (237/(1e4))*(kx*cos(th1_sol1) + ky*sin(th1_sol1));
F = D - A*sin(th3_sol1);

if (C~=0) && (B~=0)
    th2_sol1 = sign(C)*asin(F/sqrt(B^2 + C^2)) - sign(B)*sign(C)*atan(abs(B/C)) - th3_sol1;
elseif (C == 0) 
    th2_sol1 = acos(F/B) - th3_sol1;
elseif (B == 0)
    th2_sol1 = asin(F/C) - th3_sol1;
end

% CALCULATION OF THETA2 (SOLUTION 2)

B = px*cos(th1_sol1) + py*sin(th1_sol1) - (237/(1e4))*(kx*cos(th1_sol1) + ky*sin(th1_sol1));
F = D - A*sin(th3_sol1);

if (C~=0) && (B~=0)
   th2_sol2 = (pi - sign(C)*asin(F/sqrt(B^2 + C^2))) - sign(B)*sign(C)*atan(abs(B/C)) - th3_sol1;
elseif (C == 0) 
   th2_sol2 = -acos(F/B) - th3_sol1;
elseif (B == 0)
   th2_sol2 = (pi - asin(F/C)) - th3_sol1;
end

% CALCULATION OF THETA2 (SOLUTION 3)

B = px*cos(th1_sol2) + py*sin(th1_sol2) - (237/(1e4))*(kx*cos(th1_sol2) + ky*sin(th1_sol2));
F = D - A*sin(th3_sol2);

if (C~=0) && (B~=0)
    th2_sol3 = sign(C)*asin(F/sqrt(B^2 + C^2)) - sign(B)*sign(C)*atan(abs(B/C)) - th3_sol2;
elseif (C == 0) 
    th2_sol3 = acos(F/B) - th3_sol2;
elseif (B == 0)
    th2_sol3 = asin(F/C) - th3_sol2;
end

% CALCULATION OF THETA2 (SOLUTION 4)

B = px*cos(th1_sol2) + py*sin(th1_sol2) - (237/(1e4))*(kx*cos(th1_sol2) + ky*sin(th1_sol2));
F = D - A*sin(th3_sol2);

if (C~=0) && (B~=0)
   th2_sol4 = (pi - sign(C)*asin(F/sqrt(B^2 + C^2))) - sign(B)*sign(C)*atan(abs(B/C)) - th3_sol2;
elseif (C == 0) 
   th2_sol4 = -acos(F/B) - th3_sol2;
elseif (B == 0)
   th2_sol4 = (pi - asin(F/C)) - th3_sol2;
end

%% THETA5 - 8 SOLUTIONS 

% CALCULATION OF THETA5 (SOLUTION 1)

th5_sol1 = -acos(kx*cos(th1_sol1)*cos(th2_sol1)*cos(th3_sol1) - kz*cos(th3_sol1)*sin(th2_sol1) - kz*cos(th2_sol1)*sin(th3_sol1) + ky*cos(th2_sol1)*cos(th3_sol1)*sin(th1_sol1) - kx*cos(th1_sol1)*sin(th2_sol1)*sin(th3_sol1) - ky*sin(th1_sol1)*sin(th2_sol1)*sin(th3_sol1));

% CALCULATION OF THETA5 (SOLUTION 2)

th5_sol2 = acos(kx*cos(th1_sol1)*cos(th2_sol1)*cos(th3_sol1) - kz*cos(th3_sol1)*sin(th2_sol1) - kz*cos(th2_sol1)*sin(th3_sol1) + ky*cos(th2_sol1)*cos(th3_sol1)*sin(th1_sol1) - kx*cos(th1_sol1)*sin(th2_sol1)*sin(th3_sol1) - ky*sin(th1_sol1)*sin(th2_sol1)*sin(th3_sol1));

% CALCULATION OF THETA5 (SOLUTION 3)

th5_sol3 = -acos(kx*cos(th1_sol1)*cos(th2_sol2)*cos(th3_sol1) - kz*cos(th3_sol1)*sin(th2_sol2) - kz*cos(th2_sol2)*sin(th3_sol1) + ky*cos(th2_sol2)*cos(th3_sol1)*sin(th1_sol1) - kx*cos(th1_sol1)*sin(th2_sol2)*sin(th3_sol1) - ky*sin(th1_sol1)*sin(th2_sol2)*sin(th3_sol1));

% CALCULATION OF THETA5 (SOLUTION 4)

th5_sol4 = acos(kx*cos(th1_sol1)*cos(th2_sol2)*cos(th3_sol1) - kz*cos(th3_sol1)*sin(th2_sol2) - kz*cos(th2_sol2)*sin(th3_sol1) + ky*cos(th2_sol2)*cos(th3_sol1)*sin(th1_sol1) - kx*cos(th1_sol1)*sin(th2_sol2)*sin(th3_sol1) - ky*sin(th1_sol1)*sin(th2_sol2)*sin(th3_sol1));

% CALCULATION OF THETA5 (SOLUTION 5)

th5_sol5 = -acos(kx*cos(th1_sol2)*cos(th2_sol3)*cos(th3_sol2) - kz*cos(th3_sol2)*sin(th2_sol3) - kz*cos(th2_sol3)*sin(th3_sol2) + ky*cos(th2_sol3)*cos(th3_sol2)*sin(th1_sol2) - kx*cos(th1_sol2)*sin(th2_sol3)*sin(th3_sol2) - ky*sin(th1_sol2)*sin(th2_sol3)*sin(th3_sol2));

% CALCULATION OF THETA5 (SOLUTION 6)

th5_sol6 = acos(kx*cos(th1_sol2)*cos(th2_sol3)*cos(th3_sol2) - kz*cos(th3_sol2)*sin(th2_sol3) - kz*cos(th2_sol3)*sin(th3_sol2) + ky*cos(th2_sol3)*cos(th3_sol2)*sin(th1_sol2) - kx*cos(th1_sol2)*sin(th2_sol3)*sin(th3_sol2) - ky*sin(th1_sol2)*sin(th2_sol3)*sin(th3_sol2));

% CALCULATION OF THETA5 (SOLUTION 7)

th5_sol7 = -acos(kx*cos(th1_sol2)*cos(th2_sol4)*cos(th3_sol2) - kz*cos(th3_sol2)*sin(th2_sol4) - kz*cos(th2_sol4)*sin(th3_sol2) + ky*cos(th2_sol4)*cos(th3_sol2)*sin(th1_sol2) - kx*cos(th1_sol2)*sin(th2_sol4)*sin(th3_sol2) - ky*sin(th1_sol2)*sin(th2_sol4)*sin(th3_sol2));

% CALCULATION OF THETA5 (SOLUTION 8)

th5_sol8 = acos(kx*cos(th1_sol2)*cos(th2_sol4)*cos(th3_sol2) - kz*cos(th3_sol2)*sin(th2_sol4) - kz*cos(th2_sol4)*sin(th3_sol2) + ky*cos(th2_sol4)*cos(th3_sol2)*sin(th1_sol2) - kx*cos(th1_sol2)*sin(th2_sol4)*sin(th3_sol2) - ky*sin(th1_sol2)*sin(th2_sol4)*sin(th3_sol2));

%% THETA4 - 16 SOLUTIONS 

% CALCULATION OF THETA4 (SOLUTION 1)

if (th5_sol1 ~= 0)
    th4_sol1 = asin((ky*cos(th1_sol1) - kx*sin(th1_sol1))/sin(th5_sol1));
else
    th4_sol1 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 2)

if (th5_sol1 ~= 0)
    th4_sol2 = pi - asin((ky*cos(th1_sol1) - kx*sin(th1_sol1))/sin(th5_sol1));
else
    th4_sol2 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 3)

if (th5_sol2 ~= 0)
    th4_sol3 = asin((ky*cos(th1_sol1) - kx*sin(th1_sol1))/sin(th5_sol2));
else
    th4_sol3 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 4)

if (th5_sol2 ~= 0)
    th4_sol4 = pi - asin((ky*cos(th1_sol1) - kx*sin(th1_sol1))/sin(th5_sol2));
else
    th4_sol4 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 5)

if (th5_sol3 ~= 0)
    th4_sol5 = asin((ky*cos(th1_sol1) - kx*sin(th1_sol1))/sin(th5_sol3));
else
    th4_sol5 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 6)

if (th5_sol3 ~= 0)
    th4_sol6 = pi - asin((ky*cos(th1_sol1) - kx*sin(th1_sol1))/sin(th5_sol3));
else
    th4_sol6 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 7)

if (th5_sol4 ~= 0)
    th4_sol7 = asin((ky*cos(th1_sol1) - kx*sin(th1_sol1))/sin(th5_sol4));
else
    th4_sol7 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 8)

if (th5_sol4 ~= 0)
    th4_sol8 = pi - asin((ky*cos(th1_sol1) - kx*sin(th1_sol1))/sin(th5_sol4));
else
    th4_sol8 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 9)

if (th5_sol5 ~= 0)
    th4_sol9 = asin((ky*cos(th1_sol2) - kx*sin(th1_sol2))/sin(th5_sol5));
else
    th4_sol9 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 10)

if (th5_sol5 ~= 0)
    th4_sol10 = pi - asin((ky*cos(th1_sol2) - kx*sin(th1_sol2))/sin(th5_sol5));
else
    th4_sol10 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 11)

if (th5_sol6 ~= 0)
    th4_sol11 = asin((ky*cos(th1_sol2) - kx*sin(th1_sol2))/sin(th5_sol6));
else
    th4_sol11 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 12)

if (th5_sol6 ~= 0)
    th4_sol12 = pi - asin((ky*cos(th1_sol2) - kx*sin(th1_sol2))/sin(th5_sol6));
else
    th4_sol12 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 13)

if (th5_sol7 ~= 0)
    th4_sol13 = asin((ky*cos(th1_sol2) - kx*sin(th1_sol2))/sin(th5_sol7));
else
    th4_sol13 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 14)

if (th5_sol7 ~= 0)
    th4_sol14 = pi - asin((ky*cos(th1_sol2) - kx*sin(th1_sol2))/sin(th5_sol7));
else
    th4_sol14 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 15)

if (th5_sol8 ~= 0)
    th4_sol15 = asin((ky*cos(th1_sol2) - kx*sin(th1_sol2))/sin(th5_sol8));
else
    th4_sol15 = 0;
end

% CALCULATION OF THETA4 (SOLUTION 16)

if (th5_sol8 ~= 0)
    th4_sol16 = pi - asin((ky*cos(th1_sol2) - kx*sin(th1_sol2))/sin(th5_sol8));
else
    th4_sol16 = 0;
end

%% THETA6 - 32 SOLUTIONS 

% CALCULATION OF THETA6 (SOLUTION 1)

th6_sol1 = asin(ix*cos(th4_sol1)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol1) - iz*cos(th2_sol1)*cos(th3_sol1)*sin(th4_sol1) + iz*sin(th2_sol1)*sin(th3_sol1)*sin(th4_sol1) - ix*cos(th1_sol1)*cos(th2_sol1)*sin(th3_sol1)*sin(th4_sol1) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol1)*sin(th4_sol1) - iy*cos(th2_sol1)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol1) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol1)*sin(th4_sol1));

% CALCULATION OF THETA6 (SOLUTION 2)

th6_sol2 = pi - asin(ix*cos(th4_sol1)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol1) - iz*cos(th2_sol1)*cos(th3_sol1)*sin(th4_sol1) + iz*sin(th2_sol1)*sin(th3_sol1)*sin(th4_sol1) - ix*cos(th1_sol1)*cos(th2_sol1)*sin(th3_sol1)*sin(th4_sol1) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol1)*sin(th4_sol1) - iy*cos(th2_sol1)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol1) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol1)*sin(th4_sol1));

% CALCULATION OF THETA6 (SOLUTION 3)

th6_sol3 = asin(ix*cos(th4_sol2)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol2) - iz*cos(th2_sol1)*cos(th3_sol1)*sin(th4_sol2) + iz*sin(th2_sol1)*sin(th3_sol1)*sin(th4_sol2) - ix*cos(th1_sol1)*cos(th2_sol1)*sin(th3_sol1)*sin(th4_sol2) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol1)*sin(th4_sol2) - iy*cos(th2_sol1)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol2) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol1)*sin(th4_sol2));

% CALCULATION OF THETA6 (SOLUTION 4)

th6_sol4 = pi - asin(ix*cos(th4_sol2)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol2) - iz*cos(th2_sol1)*cos(th3_sol1)*sin(th4_sol2) + iz*sin(th2_sol1)*sin(th3_sol1)*sin(th4_sol2) - ix*cos(th1_sol1)*cos(th2_sol1)*sin(th3_sol1)*sin(th4_sol2) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol1)*sin(th4_sol2) - iy*cos(th2_sol1)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol2) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol1)*sin(th4_sol2));

% CALCULATION OF THETA6 (SOLUTION 5)

th6_sol5 = asin(ix*cos(th4_sol3)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol3) - iz*cos(th2_sol1)*cos(th3_sol1)*sin(th4_sol3) + iz*sin(th2_sol1)*sin(th3_sol1)*sin(th4_sol3) - ix*cos(th1_sol1)*cos(th2_sol1)*sin(th3_sol1)*sin(th4_sol3) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol1)*sin(th4_sol3) - iy*cos(th2_sol1)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol3) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol1)*sin(th4_sol3));

% CALCULATION OF THETA6 (SOLUTION 6)

th6_sol6 = pi - asin(ix*cos(th4_sol3)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol3) - iz*cos(th2_sol1)*cos(th3_sol1)*sin(th4_sol3) + iz*sin(th2_sol1)*sin(th3_sol1)*sin(th4_sol3) - ix*cos(th1_sol1)*cos(th2_sol1)*sin(th3_sol1)*sin(th4_sol3) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol1)*sin(th4_sol3) - iy*cos(th2_sol1)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol3) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol1)*sin(th4_sol3));

% CALCULATION OF THETA6 (SOLUTION 7)

th6_sol7 = asin(ix*cos(th4_sol4)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol4) - iz*cos(th2_sol1)*cos(th3_sol1)*sin(th4_sol4) + iz*sin(th2_sol1)*sin(th3_sol1)*sin(th4_sol4) - ix*cos(th1_sol1)*cos(th2_sol1)*sin(th3_sol1)*sin(th4_sol4) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol1)*sin(th4_sol4) - iy*cos(th2_sol1)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol4) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol1)*sin(th4_sol4));

% CALCULATION OF THETA6 (SOLUTION 8)

th6_sol8 = pi - asin(ix*cos(th4_sol4)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol4) - iz*cos(th2_sol1)*cos(th3_sol1)*sin(th4_sol4) + iz*sin(th2_sol1)*sin(th3_sol1)*sin(th4_sol4) - ix*cos(th1_sol1)*cos(th2_sol1)*sin(th3_sol1)*sin(th4_sol4) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol1)*sin(th4_sol4) - iy*cos(th2_sol1)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol4) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol1)*sin(th4_sol4));

% CALCULATION OF THETA6 (SOLUTION 9)

th6_sol9 = asin(ix*cos(th4_sol5)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol5) - iz*cos(th2_sol2)*cos(th3_sol1)*sin(th4_sol5) + iz*sin(th2_sol2)*sin(th3_sol1)*sin(th4_sol5) - ix*cos(th1_sol1)*cos(th2_sol2)*sin(th3_sol1)*sin(th4_sol5) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol2)*sin(th4_sol5) - iy*cos(th2_sol2)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol5) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol2)*sin(th4_sol5));

% CALCULATION OF THETA6 (SOLUTION 10)

th6_sol10 = pi - asin(ix*cos(th4_sol5)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol5) - iz*cos(th2_sol2)*cos(th3_sol1)*sin(th4_sol5) + iz*sin(th2_sol2)*sin(th3_sol1)*sin(th4_sol5) - ix*cos(th1_sol1)*cos(th2_sol2)*sin(th3_sol1)*sin(th4_sol5) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol2)*sin(th4_sol5) - iy*cos(th2_sol2)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol5) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol2)*sin(th4_sol5));

% CALCULATION OF THETA6 (SOLUTION 11)

th6_sol11 = asin(ix*cos(th4_sol6)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol6) - iz*cos(th2_sol2)*cos(th3_sol1)*sin(th4_sol6) + iz*sin(th2_sol2)*sin(th3_sol1)*sin(th4_sol6) - ix*cos(th1_sol1)*cos(th2_sol2)*sin(th3_sol1)*sin(th4_sol6) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol2)*sin(th4_sol6) - iy*cos(th2_sol2)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol6) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol2)*sin(th4_sol6));

% CALCULATION OF THETA6 (SOLUTION 12)

th6_sol12 = pi - asin(ix*cos(th4_sol6)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol6) - iz*cos(th2_sol2)*cos(th3_sol1)*sin(th4_sol6) + iz*sin(th2_sol2)*sin(th3_sol1)*sin(th4_sol6) - ix*cos(th1_sol1)*cos(th2_sol2)*sin(th3_sol1)*sin(th4_sol6) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol2)*sin(th4_sol6) - iy*cos(th2_sol2)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol6) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol2)*sin(th4_sol6));

% CALCULATION OF THETA6 (SOLUTION 13)

th6_sol13 = asin(ix*cos(th4_sol7)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol7) - iz*cos(th2_sol2)*cos(th3_sol1)*sin(th4_sol7) + iz*sin(th2_sol2)*sin(th3_sol1)*sin(th4_sol7) - ix*cos(th1_sol1)*cos(th2_sol2)*sin(th3_sol1)*sin(th4_sol7) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol2)*sin(th4_sol7) - iy*cos(th2_sol2)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol7) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol2)*sin(th4_sol7));

% CALCULATION OF THETA6 (SOLUTION 14)

th6_sol14 = pi - asin(ix*cos(th4_sol7)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol7) - iz*cos(th2_sol2)*cos(th3_sol1)*sin(th4_sol7) + iz*sin(th2_sol2)*sin(th3_sol1)*sin(th4_sol7) - ix*cos(th1_sol1)*cos(th2_sol2)*sin(th3_sol1)*sin(th4_sol7) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol2)*sin(th4_sol7) - iy*cos(th2_sol2)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol7) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol2)*sin(th4_sol7));

% CALCULATION OF THETA6 (SOLUTION 15)

th6_sol15 = asin(ix*cos(th4_sol8)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol8) - iz*cos(th2_sol2)*cos(th3_sol1)*sin(th4_sol8) + iz*sin(th2_sol2)*sin(th3_sol1)*sin(th4_sol8) - ix*cos(th1_sol1)*cos(th2_sol2)*sin(th3_sol1)*sin(th4_sol8) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol2)*sin(th4_sol8) - iy*cos(th2_sol2)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol8) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol2)*sin(th4_sol8));

% CALCULATION OF THETA6 (SOLUTION 16)

th6_sol16 = pi - asin(ix*cos(th4_sol8)*sin(th1_sol1)- iy*cos(th1_sol1)*cos(th4_sol8) - iz*cos(th2_sol2)*cos(th3_sol1)*sin(th4_sol8) + iz*sin(th2_sol2)*sin(th3_sol1)*sin(th4_sol8) - ix*cos(th1_sol1)*cos(th2_sol2)*sin(th3_sol1)*sin(th4_sol8) - ix*cos(th1_sol1)*cos(th3_sol1)*sin(th2_sol2)*sin(th4_sol8) - iy*cos(th2_sol2)*sin(th1_sol1)*sin(th3_sol1)*sin(th4_sol8) - iy*cos(th3_sol1)*sin(th1_sol1)*sin(th2_sol2)*sin(th4_sol8));

% CALCULATION OF THETA6 (SOLUTION 17)

th6_sol17 = asin(ix*cos(th4_sol9)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol9) - iz*cos(th2_sol3)*cos(th3_sol2)*sin(th4_sol9) + iz*sin(th2_sol3)*sin(th3_sol2)*sin(th4_sol9) - ix*cos(th1_sol2)*cos(th2_sol3)*sin(th3_sol2)*sin(th4_sol9) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol3)*sin(th4_sol9) - iy*cos(th2_sol3)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol9) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol3)*sin(th4_sol9));

% CALCULATION OF THETA6 (SOLUTION 18)

th6_sol18 = pi - asin(ix*cos(th4_sol9)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol9) - iz*cos(th2_sol3)*cos(th3_sol2)*sin(th4_sol9) + iz*sin(th2_sol3)*sin(th3_sol2)*sin(th4_sol9) - ix*cos(th1_sol2)*cos(th2_sol3)*sin(th3_sol2)*sin(th4_sol9) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol3)*sin(th4_sol9) - iy*cos(th2_sol3)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol9) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol3)*sin(th4_sol9));

% CALCULATION OF THETA6 (SOLUTION 19)

th6_sol19 = asin(ix*cos(th4_sol10)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol10) - iz*cos(th2_sol3)*cos(th3_sol2)*sin(th4_sol10) + iz*sin(th2_sol3)*sin(th3_sol2)*sin(th4_sol10) - ix*cos(th1_sol2)*cos(th2_sol3)*sin(th3_sol2)*sin(th4_sol10) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol3)*sin(th4_sol10) - iy*cos(th2_sol3)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol10) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol3)*sin(th4_sol10));

% CALCULATION OF THETA6 (SOLUTION 20)

th6_sol20 = pi - asin(ix*cos(th4_sol10)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol10) - iz*cos(th2_sol3)*cos(th3_sol2)*sin(th4_sol10) + iz*sin(th2_sol3)*sin(th3_sol2)*sin(th4_sol10) - ix*cos(th1_sol2)*cos(th2_sol3)*sin(th3_sol2)*sin(th4_sol10) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol3)*sin(th4_sol10) - iy*cos(th2_sol3)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol10) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol3)*sin(th4_sol10));

% CALCULATION OF THETA6 (SOLUTION 21)

th6_sol21 = asin(ix*cos(th4_sol11)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol11) - iz*cos(th2_sol3)*cos(th3_sol2)*sin(th4_sol11) + iz*sin(th2_sol3)*sin(th3_sol2)*sin(th4_sol11) - ix*cos(th1_sol2)*cos(th2_sol3)*sin(th3_sol2)*sin(th4_sol11) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol3)*sin(th4_sol11) - iy*cos(th2_sol3)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol11) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol3)*sin(th4_sol11));

% CALCULATION OF THETA6 (SOLUTION 22)

th6_sol22 = pi - asin(ix*cos(th4_sol11)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol11) - iz*cos(th2_sol3)*cos(th3_sol2)*sin(th4_sol11) + iz*sin(th2_sol3)*sin(th3_sol2)*sin(th4_sol11) - ix*cos(th1_sol2)*cos(th2_sol3)*sin(th3_sol2)*sin(th4_sol11) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol3)*sin(th4_sol11) - iy*cos(th2_sol3)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol11) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol3)*sin(th4_sol11));

% CALCULATION OF THETA6 (SOLUTION 23)

th6_sol23 = asin(ix*cos(th4_sol12)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol12) - iz*cos(th2_sol3)*cos(th3_sol2)*sin(th4_sol12) + iz*sin(th2_sol3)*sin(th3_sol2)*sin(th4_sol12) - ix*cos(th1_sol2)*cos(th2_sol3)*sin(th3_sol2)*sin(th4_sol12) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol3)*sin(th4_sol12) - iy*cos(th2_sol3)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol12) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol3)*sin(th4_sol12));

% CALCULATION OF THETA6 (SOLUTION 24)

th6_sol24 = pi - asin(ix*cos(th4_sol12)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol12) - iz*cos(th2_sol3)*cos(th3_sol2)*sin(th4_sol12) + iz*sin(th2_sol3)*sin(th3_sol2)*sin(th4_sol12) - ix*cos(th1_sol2)*cos(th2_sol3)*sin(th3_sol2)*sin(th4_sol12) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol3)*sin(th4_sol12) - iy*cos(th2_sol3)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol12) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol3)*sin(th4_sol12));

% CALCULATION OF THETA6 (SOLUTION 25)

th6_sol25 = asin(ix*cos(th4_sol13)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol13) - iz*cos(th2_sol4)*cos(th3_sol2)*sin(th4_sol13) + iz*sin(th2_sol4)*sin(th3_sol2)*sin(th4_sol13) - ix*cos(th1_sol2)*cos(th2_sol4)*sin(th3_sol2)*sin(th4_sol13) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol4)*sin(th4_sol13) - iy*cos(th2_sol4)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol13) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol4)*sin(th4_sol13));

% CALCULATION OF THETA6 (SOLUTION 26)

th6_sol26 = pi - asin(ix*cos(th4_sol13)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol13) - iz*cos(th2_sol4)*cos(th3_sol2)*sin(th4_sol13) + iz*sin(th2_sol4)*sin(th3_sol2)*sin(th4_sol13) - ix*cos(th1_sol2)*cos(th2_sol4)*sin(th3_sol2)*sin(th4_sol13) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol4)*sin(th4_sol13) - iy*cos(th2_sol4)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol13) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol4)*sin(th4_sol13));

% CALCULATION OF THETA6 (SOLUTION 27)

th6_sol27 = asin(ix*cos(th4_sol14)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol14) - iz*cos(th2_sol4)*cos(th3_sol2)*sin(th4_sol14) + iz*sin(th2_sol4)*sin(th3_sol2)*sin(th4_sol14) - ix*cos(th1_sol2)*cos(th2_sol4)*sin(th3_sol2)*sin(th4_sol14) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol4)*sin(th4_sol14) - iy*cos(th2_sol4)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol14) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol4)*sin(th4_sol14));

% CALCULATION OF THETA6 (SOLUTION 28)

th6_sol28 = pi - asin(ix*cos(th4_sol14)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol14) - iz*cos(th2_sol4)*cos(th3_sol2)*sin(th4_sol14) + iz*sin(th2_sol4)*sin(th3_sol2)*sin(th4_sol14) - ix*cos(th1_sol2)*cos(th2_sol4)*sin(th3_sol2)*sin(th4_sol14) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol4)*sin(th4_sol14) - iy*cos(th2_sol4)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol14) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol4)*sin(th4_sol14));

% CALCULATION OF THETA6 (SOLUTION 29)

th6_sol29 = asin(ix*cos(th4_sol15)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol15) - iz*cos(th2_sol4)*cos(th3_sol2)*sin(th4_sol15) + iz*sin(th2_sol4)*sin(th3_sol2)*sin(th4_sol15) - ix*cos(th1_sol2)*cos(th2_sol4)*sin(th3_sol2)*sin(th4_sol15) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol4)*sin(th4_sol15) - iy*cos(th2_sol4)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol15) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol4)*sin(th4_sol15));

% CALCULATION OF THETA6 (SOLUTION 30)

th6_sol30 = pi - asin(ix*cos(th4_sol15)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol15) - iz*cos(th2_sol4)*cos(th3_sol2)*sin(th4_sol15) + iz*sin(th2_sol4)*sin(th3_sol2)*sin(th4_sol15) - ix*cos(th1_sol2)*cos(th2_sol4)*sin(th3_sol2)*sin(th4_sol15) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol4)*sin(th4_sol15) - iy*cos(th2_sol4)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol15) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol4)*sin(th4_sol15));

% CALCULATION OF THETA6 (SOLUTION 31)

th6_sol31 = asin(ix*cos(th4_sol16)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol16) - iz*cos(th2_sol4)*cos(th3_sol2)*sin(th4_sol16) + iz*sin(th2_sol4)*sin(th3_sol2)*sin(th4_sol16) - ix*cos(th1_sol2)*cos(th2_sol4)*sin(th3_sol2)*sin(th4_sol16) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol4)*sin(th4_sol16) - iy*cos(th2_sol4)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol16) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol4)*sin(th4_sol16));

% CALCULATION OF THETA6 (SOLUTION 32)

th6_sol32 = pi - asin(ix*cos(th4_sol16)*sin(th1_sol2)- iy*cos(th1_sol2)*cos(th4_sol16) - iz*cos(th2_sol4)*cos(th3_sol2)*sin(th4_sol16) + iz*sin(th2_sol4)*sin(th3_sol2)*sin(th4_sol16) - ix*cos(th1_sol2)*cos(th2_sol4)*sin(th3_sol2)*sin(th4_sol16) - ix*cos(th1_sol2)*cos(th3_sol2)*sin(th2_sol4)*sin(th4_sol16) - iy*cos(th2_sol4)*sin(th1_sol2)*sin(th3_sol2)*sin(th4_sol16) - iy*cos(th3_sol2)*sin(th1_sol2)*sin(th2_sol4)*sin(th4_sol16));


%% REWRITE THE SOLUTIONS IN VECTOR FORM

th1 = [th1_sol1 th1_sol2];

th2 = [th2_sol1 th2_sol2 th2_sol3 th2_sol4];

th3 = [th3_sol1 th3_sol1];

th4 = [th4_sol1 th4_sol2 th4_sol3 th4_sol4 th4_sol5 th4_sol6 th4_sol7 th4_sol8 th4_sol9 th4_sol10 th4_sol11 th4_sol12 th4_sol13 th4_sol14 th4_sol15 th4_sol16];

th5 = [th5_sol1 th5_sol2 th5_sol3 th5_sol4 th5_sol5 th5_sol6 th5_sol7 th5_sol8];

th6 = [th6_sol1 th6_sol2 th6_sol3 th6_sol4 th6_sol5 th6_sol6 th6_sol7 th6_sol8 th6_sol9 th6_sol10 th6_sol11 th6_sol12 th6_sol13 th6_sol14 th6_sol15 th6_sol16 th6_sol17 th6_sol18 th6_sol19 th6_sol20 th6_sol21 th6_sol22 th6_sol23 th6_sol24 th6_sol25 th6_sol26 th6_sol27 th6_sol28 th6_sol29 th6_sol30 th6_sol31 th6_sol32];

% SOLUTION SETS (32)

THETA = zeros(32,6);

for i = 1:32
  if (i<=2)
    THETA(i,:) = [th1(1) th2(1) th3(1) th4(1) th5(1) th6(i)];
    THETA(i,:) = [th1(1) th2(1) th3(1) th4(1) th5(1) th6(i)];
  elseif (i>2) && (i<=4)
    THETA(i,:) = [th1(1) th2(1) th3(1) th4(2) th5(1) th6(i)];
    THETA(i,:) = [th1(1) th2(1) th3(1) th4(2) th5(1) th6(i)];
  elseif (i>4) && (i<=6)
    THETA(i,:) = [th1(1) th2(1) th3(1) th4(3) th5(2) th6(i)];
    THETA(i,:) = [th1(1) th2(1) th3(1) th4(3) th5(2) th6(i)];
  elseif (i>6) && (i<=8) 
    THETA(i,:) = [th1(1) th2(1) th3(1) th4(4) th5(2) th6(i)];
    THETA(i,:) = [th1(1) th2(1) th3(1) th4(4) th5(2) th6(i)];
  elseif (i>8) && (i<=10)
    THETA(i,:) = [th1(1) th2(2) th3(1) th4(5) th5(3) th6(i)];
    THETA(i,:) = [th1(1) th2(2) th3(1) th4(5) th5(3) th6(i)];
  elseif (i>10) && (i<=12)
    THETA(i,:) = [th1(1) th2(2) th3(1) th4(6) th5(3) th6(i)];
    THETA(i,:) = [th1(1) th2(2) th3(1) th4(6) th5(3) th6(i)];
  elseif (i>12) && (i<=14)
    THETA(i,:) = [th1(1) th2(2) th3(1) th4(7) th5(4) th6(i)];
    THETA(i,:) = [th1(1) th2(2) th3(1) th4(7) th5(4) th6(i)];
  elseif (i>14) && (i<=16)
    THETA(i,:) = [th1(1) th2(2) th3(1) th4(8) th5(4) th6(i)];
    THETA(i,:) = [th1(1) th2(2) th3(1) th4(8) th5(4) th6(i)];
  elseif (i>16) && (i<=18)
    THETA(i,:) = [th1(2) th2(3) th3(2) th4(9) th5(5) th6(i)];
    THETA(i,:) = [th1(2) th2(3) th3(2) th4(9) th5(5) th6(i)];
  elseif (i>18) && (i<=20)
    THETA(i,:) = [th1(2) th2(3) th3(2) th4(10) th5(5) th6(i)];
    THETA(i,:) = [th1(2) th2(3) th3(2) th4(10) th5(5) th6(i)];
  elseif (i>20) && (i<=22)
    THETA(i,:) = [th1(2) th2(3) th3(2) th4(11) th5(6) th6(i)];
    THETA(i,:) = [th1(2) th2(3) th3(2) th4(11) th5(6) th6(i)];
  elseif (i>22) && (i<=24)
    THETA(i,:) = [th1(2) th2(3) th3(2) th4(12) th5(6) th6(i)];
    THETA(i,:) = [th1(2) th2(3) th3(2) th4(12) th5(6) th6(i)];
  elseif (i>24) && (i<=26)
    THETA(i,:) = [th1(2) th2(4) th3(2) th4(13) th5(7) th6(i)];
    THETA(i,:) = [th1(2) th2(4) th3(2) th4(13) th5(7) th6(i)];
  elseif (i>26) && (i<=28)  
    THETA(i,:) = [th1(2) th2(4) th3(2) th4(14) th5(7) th6(i)];
    THETA(i,:) = [th1(2) th2(4) th3(2) th4(14) th5(7) th6(i)];
  elseif (i>28) && (i<=30)
    THETA(i,:) = [th1(2) th2(4) th3(2) th4(15) th5(8) th6(i)];
    THETA(i,:) = [th1(2) th2(4) th3(2) th4(15) th5(8) th6(i)];
  elseif (i>30) && (i<=32)
    THETA(i,:) = [th1(2) th2(4) th3(2) th4(16) th5(8) th6(i)];
    THETA(i,:) = [th1(2) th2(4) th3(2) th4(16) th5(8) th6(i)];
  end      
end

%% CHECK THE SOLUTION 

a = ([0 -90 0 -90 90 -90])*pi/180;
L = ([0 0 210 30 0 0])*(1e-3);
d = ([80 + 103 0 0 180 + 41.5 0 23.7])*(1e-3);

j = 0;
r = 1;

th_sol = zeros(1,6);

while (r==1) && (j<=31)
    j = j + 1;
    th = [THETA(j,1) (THETA(j,2)-pi/2) THETA(j,3) THETA(j,4) THETA(j,5) THETA(j,6)];

    M01 = [cos(th(1)) -sin(th(1)) 0 L(1) ; sin(th(1))*cos(a(1)) cos(th(1))*cos(a(1)) -sin(a(1)) -sin(a(1))*d(1) ; sin(th(1))*sin(a(1)) cos(th(1))*sin(a(1)) cos(a(1)) cos(a(1))*d(1) ; 0 0 0 1];

    M12 = [cos(th(2)) -sin(th(2)) 0 L(2) ; sin(th(2))*cos(a(2)) cos(th(2))*cos(a(2)) -sin(a(2)) -sin(a(2))*d(2) ; sin(th(2))*sin(a(2)) cos(th(2))*sin(a(2)) cos(a(2)) cos(a(2))*d(2) ; 0 0 0 1];

    M23 = [cos(th(3)) -sin(th(3)) 0 L(3) ; sin(th(3))*cos(a(3)) cos(th(3))*cos(a(3)) -sin(a(3)) -sin(a(3))*d(3) ; sin(th(3))*sin(a(3)) cos(th(3))*sin(a(3)) cos(a(3)) cos(a(3))*d(3) ; 0 0 0 1];

    M34 = [cos(th(4)) -sin(th(4)) 0 L(4) ; sin(th(4))*cos(a(4)) cos(th(4))*cos(a(4)) -sin(a(4)) -sin(a(4))*d(4) ; sin(th(4))*sin(a(4)) cos(th(4))*sin(a(4)) cos(a(4)) cos(a(4))*d(4) ; 0 0 0 1];

    M45 = [cos(th(5)) -sin(th(5)) 0 L(5) ; sin(th(5))*cos(a(5)) cos(th(5))*cos(a(5)) -sin(a(5)) -sin(a(5))*d(5) ; sin(th(5))*sin(a(5)) cos(th(5))*sin(a(5)) cos(a(5)) cos(a(5))*d(5) ; 0 0 0 1];

    M56 = [cos(th(6)) -sin(th(6)) 0 L(6) ; sin(th(6))*cos(a(6)) cos(th(6))*cos(a(6)) -sin(a(6)) -sin(a(6))*d(6) ; sin(th(6))*sin(a(6)) cos(th(6))*sin(a(6)) cos(a(6)) cos(a(6))*d(6) ; 0 0 0 1];

    M06calc = M01*M12*M23*M34*M45*M56;

    M06error = abs(M06calc - M06);
    
    r = 0;
    for k = 1:4
        for l = 1:4
            if (M06error(k,l)>(1e-15)) || (abs(THETA(j,1)) > 175*pi/180) || (THETA(j,2) > 36.7*pi/180) || (THETA(j,2) < -pi/2) || (THETA(j,3) > pi/2) || (THETA(j,3) <-80*pi/180) || (abs(THETA(j,4)) > 175*pi/180) || (THETA(j,5) > 110*pi/180) || (THETA(j,5) < -100*pi/180) || (abs(THETA(j,6)) > 147.5*pi/180)
                r = 1; % THEN REPEAT THE SOLUTION 
            end
        end
    end
    
    if (r == 0)
       th_sol = THETA(j,:); % PROCEED
    end 
end

if (abs(th_sol(1)) > 175*pi/180) || (th_sol(2) > 36.7*pi/180) || (th_sol(2) < -pi/2) || (th_sol(3) > pi/2) || (th_sol(3) <-80*pi/180) || (abs(th_sol(4)) > 175*pi/180) || (th_sol(5) > 110*pi/180) || (th_sol(5) < -100*pi/180) || (abs(th_sol(6)) > 147.5*pi/180) || (r == 1)
    fprintf("\n\nSOLUTION IS NOT VALID!!\n\n");
end

% PRINT RESULTS

fprintf("\ntheta_1 = %.4f\n",round(th_sol(1)*180/pi,4));

fprintf("\ntheta_2 = %.4f\n",round(th_sol(2)*180/pi,4));

fprintf("\ntheta_3 = %.4f\n",round(th_sol(3)*180/pi,4));

fprintf("\ntheta_4 = %.4f\n",round(th_sol(4)*180/pi,4));

fprintf("\ntheta_5 = %.4f\n",round(th_sol(5)*180/pi,4));

fprintf("\ntheta_6 = %.4f\n",round(th_sol(6)*180/pi,4));

disp(M06calc);

% PLOT ROBOTIC ARM AND LOCAL COORDINATE SYSTEMS 

iloc = [20e-3 0 0 1]';
jloc = [0 20e-3 0 1]';
kloc = [0 0 20e-3 1]';

iloc1 = [30e-3 0 0 1]';
jloc1 = [0 30e-3 0 1]';
kloc1 = [0 0 30e-3 1]';

iloc2 = [25e-3 0 0 1]';
jloc2 = [0 25e-3 0 1]';
kloc2 = [0 0 25e-3 1]';

% SYSTEM 1

P1 = M01*[0 0 0 1]';

i1 = M01*iloc1;
j1 = M01*jloc1;
k1 = M01*kloc1;

% SYSTEM 2

P2 = (M01*M12)*[0 0 0 1]';

i2 = M01*M12*iloc;
j2 = M01*M12*jloc;
k2 = M01*M12*kloc;

%SYSTEM 3

P3 = (M01*M12*M23)*[0 0 0 1]';

i3 = M01*M12*M23*iloc;
j3 = M01*M12*M23*jloc;
k3 = M01*M12*M23*kloc;

% SYSTEM 4

P4 = (M01*M12*M23*M34)*[0 0 0 1]';

i4 = M01*M12*M23*M34*iloc2;
j4 = M01*M12*M23*M34*jloc2;
k4 = M01*M12*M23*M34*kloc2;

% SYSTEM 5

P5 = (M01*M12*M23*M34*M45)*[0 0 0 1]';

i5 = M01*M12*M23*M34*M45*iloc;
j5 = M01*M12*M23*M34*M45*jloc;
k5 = M01*M12*M23*M34*M45*kloc;

% SYSTEM 6

P6 = (M01*M12*M23*M34*M45*M56)*[0 0 0 1]';

i6 = M01*M12*M23*M34*M45*M56*iloc;
j6 = M01*M12*M23*M34*M45*M56*jloc;
k6 = M01*M12*M23*M34*M45*M56*kloc;

figure 

line([0 20e-3],[0 0],[0 0],"Color","blue","LineWidth",2);
text(20e-3,0,0,"x_0");
line([0 0],[0 20e-3],[0 0],"Color","green","LineWidth",2);
text(0,20e-3,0,"y_0");
line([0 0],[0 0],[0 20e-3],"Color","red","LineWidth",2);
text(0,0,20e-3,"z_0");

line([0 P1(1)],[0 P1(2)],[0 P1(3)]);

line([P1(1) i1(1)],[P1(2) i1(2)],[P1(3) i1(3)],"Color","blue","LineWidth",2);
text(i1(1),i1(2),i1(3),"x_1");
line([P1(1) j1(1)],[P1(2) j1(2)],[P1(3) j1(3)],"Color","green","LineWidth",2);
text(j1(1),j1(2),j1(3),"y_1");
line([P1(1) k1(1)],[P1(2) k1(2)],[P1(3) k1(3)],"Color","red","LineWidth",2);
text(k1(1),k1(2),k1(3),"z_1");

line([P1(1) P2(1)],[P1(2) P2(2)],[P1(3) P2(3)]);

line([P2(1) i2(1)],[P2(2) i2(2)],[P2(3) i2(3)],"Color","blue","LineWidth",2);
text(i2(1),i2(2),i2(3),"x_2");
line([P2(1) j2(1)],[P2(2) j2(2)],[P2(3) j2(3)],"Color","green","LineWidth",2);
text(j2(1),j2(2),j2(3),"y_2");
line([P2(1) k2(1)],[P2(2) k2(2)],[P2(3) k2(3)],"Color","red","LineWidth",2);
text(k2(1),k2(2),k2(3),"z_2");

line([P2(1) P3(1)],[P2(2) P3(2)],[P2(3) P3(3)]);

line([P3(1) i3(1)],[P3(2) i3(2)],[P3(3) i3(3)],"Color","blue","LineWidth",2);
text(i3(1),i3(2),i3(3),"x_3");
line([P3(1) j3(1)],[P3(2) j3(2)],[P3(3) j3(3)],"Color","green","LineWidth",2);
text(j3(1),j3(2),j3(3),"y_3");
line([P3(1) k3(1)],[P3(2) k3(2)],[P3(3) k3(3)],"Color","red","LineWidth",2);
text(k3(1),k3(2),k3(3),"z_3");

p3 = [30e-3 0 0 1]';
P34 = M01*M12*M23*p3;

line([P3(1) P34(1)],[P3(2) P34(2)],[P3(3) P34(3)]);
line([P34(1) P4(1)],[P34(2) P4(2)],[P34(3) P4(3)]);

line([P4(1) i4(1)],[P4(2) i4(2)],[P4(3) i4(3)],"Color","blue","LineWidth",2);
text(i4(1),i4(2),i4(3),"x_4");
line([P4(1) j4(1)],[P4(2) j4(2)],[P4(3) j4(3)],"Color","green","LineWidth",2);
text(j4(1),j4(2),j4(3),"y_4");
line([P4(1) k4(1)],[P4(2) k4(2)],[P4(3) k4(3)],"Color","red","LineWidth",2);
text(k4(1),k4(2),k4(3),"z_4");

line([P4(1) P5(1)],[P4(2) P5(2)],[P4(3) P5(3)]);

line([P5(1) i5(1)],[P5(2) i5(2)],[P5(3) i5(3)],"Color","blue","LineWidth",2);
text(i5(1),i5(2),i5(3),"x_5");
line([P5(1) j5(1)],[P5(2) j5(2)],[P5(3) j5(3)],"Color","green","LineWidth",2);
text(j5(1),j5(2),j5(3),"y_5");
line([P5(1) k5(1)],[P5(2) k5(2)],[P5(3) k5(3)],"Color","red","LineWidth",2);
text(k5(1),k5(2),k5(3),"z_5");

line([P5(1) P6(1)],[P5(2) P6(2)],[P5(3) P6(3)]);

line([P6(1) i6(1)],[P6(2) i6(2)],[P6(3) i6(3)],"Color","blue","LineWidth",2);
text(i6(1),i6(2),i6(3),"x_6");
line([P6(1) j6(1)],[P6(2) j6(2)],[P6(3) j6(3)],"Color","green","LineWidth",2);
text(j6(1),j6(2),j6(3),"y_6");
line([P6(1) k6(1)],[P6(2) k6(2)],[P6(3) k6(3)],"Color","red","LineWidth",2);
text(k6(1),k6(2),k6(3),"z_6");

view(3);
grid on 
xlabel("x (m)");
ylabel("y (m)");
zlabel("z (m)");
axis("equal");


