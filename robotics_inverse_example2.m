clear;
clc;

% ORIENTATION AND DISPLACEMENT PARAMETERS 

ix = 0;
iy = -1;
iz = 0;

jx = 0;
jy = 0;
jz = -1;

kx = 1;
ky = 0;
kz = 0;

px = (180 + 41.5)*(1e-3)*cos(pi/6) + 23.7e-3;
py = (180 + 41.5)*(1e-3)*sin(pi/6);
pz = 0.423;

% POSITION AND ORIENTATION OF THE WORK EDGE (SYSTEM X6Y6Z6)

M06 = [ix jx kx px ; iy jy ky py ; iz jz kz pz ; 0 0 0 1];

% CALCULATION OF THETA1 

if (-237*ky+(1e4)*py)*(-237*kx+(1e4)*px) >= 0 
    th1 = atan((-237*ky+(1e4)*py)/(-237*kx+(1e4)*px));
else 
    th1 = -atan(abs((-237*ky+(1e4)*py)/(-237*kx+(1e4)*px)));
end

% CALCULATION OF THETA3

A = 21/100;
B = px*cos(th1) + py*sin(th1) - (237/(1e4))*(kx*cos(th1) + ky*sin(th1));
C = 183/(1e3) + (237*kz/(1e4)) - pz;

D = 443/2000;
E = 3/100;

th3 = asin(((A^2 + D^2 + E^2) - (B^2 + C^2))/(2*A*sqrt(D^2 + E^2)))+atan(E/D);

% CALCULATION OF THETA2

F = D - A*sin(th3);

if C>=0
    th2 = asin(F/sqrt(B^2 + C^2)) - atan(B/C) - th3;
else
    th2 = -asin(F/sqrt(B^2 + C^2)) - atan(B/C) - th3;
end

% CALCULATION OF THETA5 (NEGATIVE SOLUTION)

th5 = -acos(kx*cos(th1)*cos(th2)*cos(th3) - kz*cos(th3)*sin(th2) - kz*cos(th2)*sin(th3) + ky*cos(th2)*cos(th3)*sin(th1) - kx*cos(th1)*sin(th2)*sin(th3) - ky*sin(th1)*sin(th2)*sin(th3));

% CALCULATION OF THETA4

if (th5 ~= 0)
    th4 = asin((ky*cos(th1) - kx*sin(th1))/sin(th5));
else
    G = (py*cos(th1) - px*sin(th1)) - (237/(1e4))*(ky*cos(th1) - kx*sin(th1));
    if (((B*sin(th2 + th3) - C*cos(th2 + th3) -A*cos(th3) -E)) == 0) && (G == 0)
        th4 = 0;
    else
        th4 = atan((B*sin(th2 + th3) - C*cos(th2 + th3) -A*cos(th3) -E)/G);
    end
end

% CALCULATION OF THETA6

th6 = asin(ix*cos(th4)*sin(th1)- iy*cos(th1)*cos(th4) - iz*cos(th2)*cos(th3)*sin(th4) + iz*sin(th2)*sin(th3)*sin(th4) - ix*cos(th1)*cos(th2)*sin(th3)*sin(th4) - ix*cos(th1)*cos(th3)*sin(th2)*sin(th4) - iy*cos(th2)*sin(th1)*sin(th3)*sin(th4) - iy*cos(th3)*sin(th1)*sin(th2)*sin(th4));

% CHECK THE SOLUTION 

a = ([0 -90 0 -90 90 -90]')*pi/180;
L = ([0 0 210 30 0 0]')*(1e-3);
th = ([th1 (th2-pi/2) th3 th4 th5 th6]');
d = ([80 + 103 0 0 180 + 41.5 0 23.7])*(1e-3);

M01 = [cos(th(1)) -sin(th(1)) 0 L(1) ; sin(th(1))*cos(a(1)) cos(th(1))*cos(a(1)) -sin(a(1)) -sin(a(1))*d(1) ; sin(th(1))*sin(a(1)) cos(th(1))*sin(a(1)) cos(a(1)) cos(a(1))*d(1) ; 0 0 0 1];

M12 = [cos(th(2)) -sin(th(2)) 0 L(2) ; sin(th(2))*cos(a(2)) cos(th(2))*cos(a(2)) -sin(a(2)) -sin(a(2))*d(2) ; sin(th(2))*sin(a(2)) cos(th(2))*sin(a(2)) cos(a(2)) cos(a(2))*d(2) ; 0 0 0 1];

M23 = [cos(th(3)) -sin(th(3)) 0 L(3) ; sin(th(3))*cos(a(3)) cos(th(3))*cos(a(3)) -sin(a(3)) -sin(a(3))*d(3) ; sin(th(3))*sin(a(3)) cos(th(3))*sin(a(3)) cos(a(3)) cos(a(3))*d(3) ; 0 0 0 1];

M34 = [cos(th(4)) -sin(th(4)) 0 L(4) ; sin(th(4))*cos(a(4)) cos(th(4))*cos(a(4)) -sin(a(4)) -sin(a(4))*d(4) ; sin(th(4))*sin(a(4)) cos(th(4))*sin(a(4)) cos(a(4)) cos(a(4))*d(4) ; 0 0 0 1];

M45 = [cos(th(5)) -sin(th(5)) 0 L(5) ; sin(th(5))*cos(a(5)) cos(th(5))*cos(a(5)) -sin(a(5)) -sin(a(5))*d(5) ; sin(th(5))*sin(a(5)) cos(th(5))*sin(a(5)) cos(a(5)) cos(a(5))*d(5) ; 0 0 0 1];

M56 = [cos(th(6)) -sin(th(6)) 0 L(6) ; sin(th(6))*cos(a(6)) cos(th(6))*cos(a(6)) -sin(a(6)) -sin(a(6))*d(6) ; sin(th(6))*sin(a(6)) cos(th(6))*sin(a(6)) cos(a(6)) cos(a(6))*d(6) ; 0 0 0 1];

M06calc = M01*M12*M23*M34*M45*M56;

M06error = abs(M06calc - M06);

r = 0;
for i = 1:4
    for j = 1:4
        if (M06error(i,j)>1e-15)
            r = 1; % THEN REPEAT THE SOLUTION WITH THE POSITIVE VALUE FOR THETA5
        end
    end
end

if (r == 1)
    th5 = acos(kx*cos(th1)*cos(th2)*cos(th3) - kz*cos(th3)*sin(th2) - kz*cos(th2)*sin(th3) + ky*cos(th2)*cos(th3)*sin(th1) - kx*cos(th1)*sin(th2)*sin(th3) - ky*sin(th1)*sin(th2)*sin(th3));

    if (th5 ~= 0)
        th4 = asin((ky*cos(th1) - kx*sin(th1))/sin(th5));
    else
        G = (py*cos(th1) - px*sin(th1)) - (237/(1e4))*(ky*cos(th1) - kx*sin(th1));
        if (((B*sin(th2 + th3) - C*cos(th2 + th3) -A*cos(th3) -E)/G) == 0) && (G == 0)
            th4 = 0;
        else
            th4 = atan((B*sin(th2 + th3) - C*cos(th2 + th3) -A*cos(th3) -E)/G);
        end
    end

    th6 = asin(ix*cos(th4)*sin(th1)- iy*cos(th1)*cos(th4) - iz*cos(th2)*cos(th3)*sin(th4) + iz*sin(th2)*sin(th3)*sin(th4) - ix*cos(th1)*cos(th2)*sin(th3)*sin(th4) - ix*cos(th1)*cos(th3)*sin(th2)*sin(th4) - iy*cos(th2)*sin(th1)*sin(th3)*sin(th4) - iy*cos(th3)*sin(th1)*sin(th2)*sin(th4));

    a = ([0 -90 0 -90 90 -90]')*pi/180;
    L = ([0 0 210 30 0 0]')*(1e-3);
    th = ([th1 (th2-pi/2) th3 th4 th5 th6]');
    d = ([80 + 103 0 0 180 + 41.5 0 23.7])*(1e-3);

    % transform matrices 
    M01 = [cos(th(1)) -sin(th(1)) 0 L(1) ; sin(th(1))*cos(a(1)) cos(th(1))*cos(a(1)) -sin(a(1)) -sin(a(1))*d(1) ; sin(th(1))*sin(a(1)) cos(th(1))*sin(a(1)) cos(a(1)) cos(a(1))*d(1) ; 0 0 0 1];
    
    M12 = [cos(th(2)) -sin(th(2)) 0 L(2) ; sin(th(2))*cos(a(2)) cos(th(2))*cos(a(2)) -sin(a(2)) -sin(a(2))*d(2) ; sin(th(2))*sin(a(2)) cos(th(2))*sin(a(2)) cos(a(2)) cos(a(2))*d(2) ; 0 0 0 1];

    M23 = [cos(th(3)) -sin(th(3)) 0 L(3) ; sin(th(3))*cos(a(3)) cos(th(3))*cos(a(3)) -sin(a(3)) -sin(a(3))*d(3) ; sin(th(3))*sin(a(3)) cos(th(3))*sin(a(3)) cos(a(3)) cos(a(3))*d(3) ; 0 0 0 1];

    M34 = [cos(th(4)) -sin(th(4)) 0 L(4) ; sin(th(4))*cos(a(4)) cos(th(4))*cos(a(4)) -sin(a(4)) -sin(a(4))*d(4) ; sin(th(4))*sin(a(4)) cos(th(4))*sin(a(4)) cos(a(4)) cos(a(4))*d(4) ; 0 0 0 1];

    M45 = [cos(th(5)) -sin(th(5)) 0 L(5) ; sin(th(5))*cos(a(5)) cos(th(5))*cos(a(5)) -sin(a(5)) -sin(a(5))*d(5) ; sin(th(5))*sin(a(5)) cos(th(5))*sin(a(5)) cos(a(5)) cos(a(5))*d(5) ; 0 0 0 1];

    M56 = [cos(th(6)) -sin(th(6)) 0 L(6) ; sin(th(6))*cos(a(6)) cos(th(6))*cos(a(6)) -sin(a(6)) -sin(a(6))*d(6) ; sin(th(6))*sin(a(6)) cos(th(6))*sin(a(6)) cos(a(6)) cos(a(6))*d(6) ; 0 0 0 1];

    M06calc = M01*M12*M23*M34*M45*M56;
  
end

% PRINT RESULTS

fprintf("\ntheta_1 = %.4f\n",round(th1*180/pi,4));

fprintf("\ntheta_2 = %.4f\n",round(th2*180/pi,4));

fprintf("\ntheta_3 = %.4f\n",round(th3*180/pi,4));

fprintf("\ntheta_4 = %.4f\n",round(th4*180/pi,4));

fprintf("\ntheta_5 = %.4f\n",round(th5*180/pi,4));

fprintf("\ntheta_6 = %.4f\n",round(th6*180/pi,4));

disp(M06calc);
