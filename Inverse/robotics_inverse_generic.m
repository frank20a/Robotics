clear;
clc;

% ORIENTATION AND DISPLACEMENT PARAMETERS 

ix = 0;
iy = 0;
iz = 1;

jx = 0;
jy = -1;
jz = 0;

kx = 1;
ky = 0;
kz = 0;

px = 0.2452;
py = -0.2;
pz = 0.1;

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

% CALCULATION OF THETA2 (SOLUTION 1)

F = D - A*sin(th3);

if (C~=0) && (B~=0)
    th2 = sign(C)*asin(F/sqrt(B^2 + C^2)) - sign(B)*sign(C)*atan(abs(B/C)) - th3;
elseif (C == 0) 
    th2 = acos(F/B) - th3;
elseif (B == 0)
    th2 = asin(F/C) - th3;
end

% CALCULATION OF THETA5 (SOLUTION 1)

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
            r = 1; % THEN REPEAT THE SOLUTION 
        end
    end
end

if (r == 1)
    % CALCULATION OF THETA2 (SOLUTION 2)

    F = D - A*sin(th3);

    if (C~=0) && (B~=0)
        th2 = (pi-sign(C)*asin(F/sqrt(B^2 + C^2))) - sign(B)*sign(C)*atan(abs(B/C)) - th3;
    elseif (C == 0) 
        th2 = -acos(F/B) - th3;
    elseif (B == 0)
        th2 = (pi - asin(F/C)) - th3;
    end

    % CALCULATION OF THETA5 (SOLUTION 1)

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

end

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
            r = 1; % THEN REPEAT THE SOLUTION 
        end
    end
end

if (r == 1)
    % CALCULATION OF THETA2 (SOLUTION 1)

    F = D - A*sin(th3);

    if (C~=0) && (B~=0)
        th2 = sign(C)*asin(F/sqrt(B^2 + C^2)) - sign(B)*sign(C)*atan(abs(B/C)) - th3;
    elseif (C == 0) 
        th2 = acos(F/B) - th3;
    elseif (B == 0)
        th2 = asin(F/C) - th3;
    end

    % CALCULATION OF THETA5 (SOLUTION 2)

    th5 = acos(kx*cos(th1)*cos(th2)*cos(th3) - kz*cos(th3)*sin(th2) - kz*cos(th2)*sin(th3) + ky*cos(th2)*cos(th3)*sin(th1) - kx*cos(th1)*sin(th2)*sin(th3) - ky*sin(th1)*sin(th2)*sin(th3));

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

end

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
            r = 1; % THEN REPEAT THE SOLUTION 
        end
    end
end

if (r == 1)
    % CALCULATION OF THETA2 (SOLUTION 2)

    F = D - A*sin(th3);

    if (C~=0) && (B~=0)
        th2 = (pi - sign(C)*asin(F/sqrt(B^2 + C^2))) - sign(B)*sign(C)*atan(abs(B/C)) - th3;
    elseif (C == 0) 
        th2 = -acos(F/B) - th3;
    elseif (B == 0)
        th2 = (pi - asin(F/C)) - th3;
    end

    % CALCULATION OF THETA5 (SOLUTION 2)

    th5 = acos(kx*cos(th1)*cos(th2)*cos(th3) - kz*cos(th3)*sin(th2) - kz*cos(th2)*sin(th3) + ky*cos(th2)*cos(th3)*sin(th1) - kx*cos(th1)*sin(th2)*sin(th3) - ky*sin(th1)*sin(th2)*sin(th3));

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

end

% FINAL CHECK 

r = 0;
for i = 1:4
    for j = 1:4
        if (M06error(i,j)>1e-15)
            r = 1; % THEN REPEAT THE SOLUTION 
        end
    end
end

if (abs(th1) > pi/2) || (abs(th2) > pi/2) || (abs(th3) > pi/2) || (abs(th4) > pi/2) || (abs(th5) > pi/2) || (abs(th6) > pi/2) || (r == 1)
    fprintf("\n\nSOLUTION IS NOT VALID!!\n\n");
end

% PRINT RESULTS

fprintf("\ntheta_1 = %.4f\n",round(th1*180/pi,4));

fprintf("\ntheta_2 = %.4f\n",round(th2*180/pi,4));

fprintf("\ntheta_3 = %.4f\n",round(th3*180/pi,4));

fprintf("\ntheta_4 = %.4f\n",round(th4*180/pi,4));

fprintf("\ntheta_5 = %.4f\n",round(th5*180/pi,4));

fprintf("\ntheta_6 = %.4f\n",round(th6*180/pi,4));

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


