clear;
clc;

%% DATA INPUT 

% FINAL POINT 
px4 = (170 - 83.62)*1e-3;
py4 = 220e-3;
pz4 = 300e-3;

% FIRST POINT 
px1 = (287.5-83.62)*(1e-3);
py1 = 62.5e-3;
pz1 = 275e-3;

% 1ST INTERMEDIATE POINT 
px2 = px1;
py2 = py1;
pz2 = pz4;

% 2ND INTERMEDIATE POINT 
px3 = px4;
py3 = py1;
pz3 = pz4;

% FIRST AND 2ND POINT ORIENTATION 
ix1 = 0;
iy1 = 0;
iz1 = 1;

jx1 = 0;
jy1 = -1;
jz1 = 0;

kx1 = 1;
ky1 = 0;
kz1 = 0;

%3RD AND FINAL POINT ORIENTATION 
ix3 = 0; 
iy3 = -1;
iz3 = 0;

jx3 = 0;
jy3 = 0;
jz3 = -1;

kx3 = 1;
ky3 = 0;
kz3 = 0;

% TRANSFORMS TO BASE 
M1 = [ix1 jx1 kx1 px1 ; iy1 jy1 ky1 py1 ; iz1 jz1 kz1 pz1 ; 0 0 0 1];

M2 = [ix1 jx1 kx1 px2 ; iy1 jy1 ky1 py2 ; iz1 jz1 kz1 pz2 ; 0 0 0 1];

M3 = [ix3 jx3 kx3 px3 ; iy3 jy3 ky3 py3 ; iz3 jz3 kz3 pz3 ; 0 0 0 1];

M4 = [ix3 jx3 kx3 px4 ; iy3 jy3 ky3 py4 ; iz3 jz3 kz3 pz4 ; 0 0 0 1];

% TIME FROM FIRST TO 2ND
t2 = 2;
% TIME FROM 2ND TO 3RD
t3 = 4;
% TIME FROM 3RD TO FINAL
t4 = 6;

% POLYNOMIAL COEFFICIENTS 
[A,B,C,D] = orbit(M1,M2,M3,M4,t2,t3,t4);

n = 1e3;

% TIME ARRAY
T = linspace(0,t4,n);

% KARTESIAN COORDINATES OF WORKING EDGE

x = zeros(1,n);
y = zeros(1,n);
z = zeros(1,n);

% ANGLE MATRIX
TH = zeros(6,n);

% TRANSFORMS MATRIX
M1 = zeros(4,4,6,n);
M2 = zeros(4,4,6,n);
M3 = zeros(4,4,6,n);


% ROBOT D-H PARAMETERS
a = ([0 -90 0 -90 90 -90])*pi/180;
L = ([0 0 210 30 0 0])*(1e-3);
d = ([80 + 103 0 0 180 + 41.5 0 23.7])*(1e-3);

% CALCULATION OF WORKING EDGE ORBIT 
for i = 1:n
    % 1ST POLYNOMIAL 
    if (T(i) <= t2)
        for j = 1:6
            TH(j,i) = A(1,j) + B(1,j)*T(i) + C(1,j)*T(i)^2 + D(1,j)*T(i)^3;
        end
        
        th = [TH(1,i) (TH(2,i)-pi/2) TH(3,i) TH(4,i) TH(5,i) TH(6,i)];
        
        % TRANSFORMS 
        for j = 1:6
            M1(:,:,j,i) = [cos(th(j)) -sin(th(j)) 0 L(j) ; sin(th(j))*cos(a(j)) cos(th(j))*cos(a(j)) -sin(a(j)) -sin(a(j))*d(j) ; sin(th(j))*sin(a(j)) cos(th(j))*sin(a(j)) cos(a(j)) cos(a(j))*d(j) ; 0 0 0 1];
        end
        
        % TRANSFORM TO BASE 
        M06calc1 = M1(:,:,1,i)*M1(:,:,2,i)*M1(:,:,3,i)*M1(:,:,4,i)*M1(:,:,5,i)*M1(:,:,6,i);
        
        % WORKNG EDGE COORDINATES
        P1 = M06calc1*[0 0 0 1]';
        
        x(i) = P1(1);
        y(i) = P1(2);
        z(i) = P1(3);
        
    elseif (T(i) <= t3) && (T(i) > t2)
        %2ND POLYNOMIAL 
       for j = 1:6
            TH(j,i) = A(2,j) + B(2,j)*T(i) + C(2,j)*T(i)^2 + D(2,j)*T(i)^3;
        end
        
        th = [TH(1,i) (TH(2,i)-pi/2) TH(3,i) TH(4,i) TH(5,i) TH(6,i)];
        
        %TRANSFORMS 
        for j = 1:6
            M2(:,:,j,i) = [cos(th(j)) -sin(th(j)) 0 L(j) ; sin(th(j))*cos(a(j)) cos(th(j))*cos(a(j)) -sin(a(j)) -sin(a(j))*d(j) ; sin(th(j))*sin(a(j)) cos(th(j))*sin(a(j)) cos(a(j)) cos(a(j))*d(j) ; 0 0 0 1];
        end
        
        % TRANSFORM TO BASE 
        M06calc2 = M2(:,:,1,i)*M2(:,:,2,i)*M2(:,:,3,i)*M2(:,:,4,i)*M2(:,:,5,i)*M2(:,:,6,i);
        
        % WORKING EDGE COORDINATES 
        P2 = M06calc2*[0 0 0 1]';
        
        x(i) = P2(1);
        y(i) = P2(2);
        z(i) = P2(3);
    elseif (T(i) <= t4) && (T(i) > t3)
        %2ND POLYNOMIAL 
       for j = 1:6
            TH(j,i) = A(3,j) + B(3,j)*T(i) + C(3,j)*T(i)^2 + D(3,j)*T(i)^3;
        end
        
        th = [TH(1,i) (TH(2,i)-pi/2) TH(3,i) TH(4,i) TH(5,i) TH(6,i)];
        
        %TRANSFORMS 
        for j = 1:6
            M3(:,:,j,i) = [cos(th(j)) -sin(th(j)) 0 L(j) ; sin(th(j))*cos(a(j)) cos(th(j))*cos(a(j)) -sin(a(j)) -sin(a(j))*d(j) ; sin(th(j))*sin(a(j)) cos(th(j))*sin(a(j)) cos(a(j)) cos(a(j))*d(j) ; 0 0 0 1];
        end
        
        % TRANSFORM TO BASE 
        M06calc3 = M3(:,:,1,i)*M3(:,:,2,i)*M3(:,:,3,i)*M3(:,:,4,i)*M3(:,:,5,i)*M3(:,:,6,i);
        
        % WORKING EDGE COORDINATES 
        P3 = M06calc3*[0 0 0 1]';
        
        x(i) = P3(1);
        y(i) = P3(2);
        z(i) = P3(3);
    end
end


% PLOT ORBIT
figure
plot3(x,y,z,"LineWidth",1.5);
hold on 
plot3(px1,py1,pz1,'.',"MarkerSize",20);
hold on 
plot3(px2,py2,pz2,'.',"MarkerSize",20);
hold on 
plot3(px3,py3,pz3,'.',"MarkerSize",20);
hold on 
plot3(px4,py4,pz4,'.',"MarkerSize",20);
title("Working Edge Orbit");
xlabel("x (m)");
ylabel("y (m)");
zlabel("z (m)");
grid on 


% CALCULATION OF POLYNOMIAL COEFFICIENTS 
function [a,b,c,d] = orbit(M1,M2,M3,M4,t2,t3,t4)
    
    % INITIALIZATION 
    a = zeros(3,6);
    b = zeros(3,6);
    c = zeros(3,6);
    d = zeros(3,6);

    % SOLUTION OF THE INVERSE KINEMATIC 
    TH1 = inverse_kinematic(M1(1,1),M1(1,2),M1(1,3),M1(1,4),M1(2,1),M1(2,2),M1(2,3),M1(2,4),M1(3,1),M1(3,2),M1(3,3),M1(3,4));
    
    TH2 = inverse_kinematic(M2(1,1),M2(1,2),M2(1,3),M2(1,4),M2(2,1),M2(2,2),M2(2,3),M2(2,4),M2(3,1),M2(3,2),M2(3,3),M2(3,4));
    
    TH3 = inverse_kinematic(M3(1,1),M3(1,2),M3(1,3),M3(1,4),M3(2,1),M3(2,2),M3(2,3),M3(2,4),M3(3,1),M3(3,2),M3(3,3),M3(3,4));
    
    TH4 = inverse_kinematic(M4(1,1),M4(1,2),M4(1,3),M4(1,4),M4(2,1),M4(2,2),M4(2,3),M4(2,4),M4(3,1),M4(3,2),M4(3,3),M4(3,4));
    
    % THE FUNCTION inverse_kinematic takes the following value as input:
    % ix,jx,kx,px,iy,jy,ky,py,iz,jz,kz,pz
    
    % CALCULATION OF POLYNOMIAL COEFFICIENTS 
    a(1,:) = TH1;
    
    c(1,:) = (3*(TH3*t2^4 - TH4*t2^4 + TH1*t2^2*t3^2 - TH1*t2^2*t4^2 + 2*TH1*t3^2*t4^2 - 2*TH2*t3^2*t4^2 + TH3*t2^2*t4^2 - TH4*t2^2*t3^2 + TH1*t2*t3^3 - 2*TH1*t2^3*t3 + 2*TH1*t2^3*t4 - TH2*t2*t3^3 - 2*TH1*t3^3*t4 + 2*TH2*t3^3*t4 - 2*TH3*t2^3*t4 + 2*TH4*t2^3*t3 - TH1*t2*t3*t4^2 + TH2*t2*t3*t4^2))/(t2^2*(t2 - t3)*(t3 - t4)*(3*t2*t3 + t2*t4 - 4*t3*t4));
    
    d(1,:) =  - (3*TH3*t2^4 - 3*TH4*t2^4 + 6*TH1*t2^2*t3^2 - 4*TH1*t2^2*t4^2 - 3*TH2*t2^2*t3^2 + 2*TH1*t3^2*t4^2 + TH2*t2^2*t4^2 - 2*TH2*t3^2*t4^2 + 3*TH3*t2^2*t4^2 - 3*TH4*t2^2*t3^2 - 6*TH1*t2^3*t3 + 6*TH1*t2^3*t4 - 2*TH1*t3^3*t4 + 2*TH2*t3^3*t4 - 6*TH3*t2^3*t4 + 6*TH4*t2^3*t3 + 2*TH1*t2*t3*t4^2 - 2*TH1*t2*t3^2*t4 - 2*TH1*t2^2*t3*t4 - 2*TH2*t2*t3*t4^2 + 2*TH2*t2*t3^2*t4 + 2*TH2*t2^2*t3*t4)/(t2^3*(t2 - t3)*(t3 - t4)*(3*t2*t3 + t2*t4 - 4*t3*t4));
    
    a(2,:) =  - (3*TH1*t2^2*t3^3 - 6*TH1*t3^3*t4^2 + 2*TH2*t3^3*t4^2 + TH3*t2^3*t4^2 + 3*TH4*t2^2*t3^3 - 3*TH4*t2^3*t3^2 - 3*TH1*t2*t3^4 + 6*TH1*t3^4*t4 - 2*TH2*t3^4*t4 - 6*TH1*t2*t3^3*t4 + 2*TH3*t2^3*t3*t4 + 9*TH1*t2*t3^2*t4^2 - 3*TH1*t2^2*t3*t4^2 - 3*TH3*t2^2*t3*t4^2)/((t2 - t3)^2*(t3 - t4)*(3*t2*t3 + t2*t4 - 4*t3*t4));
    
    b(2,:) =  - (3*(3*TH1*t2^2*t3^3 - 3*TH1*t2^3*t3^2 + TH1*t2^3*t4^2 + 2*TH1*t3^3*t4^2 - 2*TH2*t3^3*t4^2 - TH3*t2^3*t4^2 - 3*TH4*t2^2*t3^3 + 3*TH4*t2^3*t3^2 - 2*TH1*t3^4*t4 + 2*TH2*t3^4*t4 + 2*TH1*t2^3*t3*t4 - 2*TH3*t2^3*t3*t4 - 3*TH1*t2^2*t3*t4^2 + 3*TH3*t2^2*t3*t4^2))/(t2*(t2 - t3)^2*(t3 - t4)*(3*t2*t3 + t2*t4 - 4*t3*t4));
    
    c(2,:) =  - (3*(TH1*t3^4 - TH2*t3^4 - TH3*t2^4 + TH4*t2^4 - 3*TH1*t3^2*t4^2 + 3*TH2*t3^2*t4^2 - 3*TH1*t2*t3^3 + 2*TH1*t2^3*t3 - 2*TH1*t2^3*t4 + TH2*t2*t3^3 + 2*TH1*t3^3*t4 + TH3*t2^3*t3 - 2*TH2*t3^3*t4 + 2*TH3*t2^3*t4 + 2*TH4*t2*t3^3 - 3*TH4*t2^3*t3 + 3*TH1*t2*t3*t4^2 - TH2*t2*t3*t4^2 - 2*TH3*t2*t3*t4^2))/(t2*(t2 - t3)^2*(t3 - t4)*(3*t2*t3 + t2*t4 - 4*t3*t4));
    
    d(2,:) = (3*TH1*t3^3 - 3*TH2*t3^3 - 3*TH3*t2^3 + 3*TH4*t2^3 - 9*TH1*t2*t3^2 + 6*TH1*t2^2*t3 + 3*TH1*t2*t4^2 - 6*TH1*t2^2*t4 + 3*TH2*t2*t3^2 - 3*TH1*t3*t4^2 - TH2*t2*t4^2 + 3*TH3*t2^2*t3 + 3*TH2*t3*t4^2 - 2*TH3*t2*t4^2 + 6*TH3*t2^2*t4 + 6*TH4*t2*t3^2 - 9*TH4*t2^2*t3 + 6*TH1*t2*t3*t4 - 2*TH2*t2*t3*t4 - 4*TH3*t2*t3*t4)/(t2*(t2 - t3)^2*(t3 - t4)*(3*t2*t3 + t2*t4 - 4*t3*t4));
    
    a(3,:) = (3*TH1*t3^3*t4^4 - 6*TH1*t3^4*t4^3 + 3*TH1*t3^5*t4^2 - 3*TH2*t3^3*t4^4 + 6*TH2*t3^4*t4^3 - 3*TH2*t3^5*t4^2 - TH3*t2^3*t4^4 - 3*TH4*t2^2*t3^5 + 3*TH4*t2^3*t3^4 - 6*TH1*t2^2*t3^2*t4^3 + 3*TH1*t2^2*t3^3*t4^2 - 6*TH3*t2^2*t3^2*t4^3 + 3*TH3*t2^3*t3^2*t4^2 + 3*TH4*t2^2*t3^3*t4^2 + 3*TH4*t2^3*t3^2*t4^2 + 4*TH4*t2*t3^5*t4 - 6*TH1*t2*t3^2*t4^4 + 12*TH1*t2*t3^3*t4^3 - 6*TH1*t2*t3^4*t4^2 + 3*TH1*t2^2*t3*t4^4 + 2*TH3*t2*t3^2*t4^4 + 2*TH3*t2^2*t3*t4^4 - 6*TH4*t2*t3^4*t4^2 + 4*TH4*t2^2*t3^4*t4 - 8*TH4*t2^3*t3^3*t4)/(t2*(t2 - t3)*(t3 - t4)^3*(3*t2*t3 + t2*t4 - 4*t3*t4));
    
    b(3,:) = (3*t4*(2*TH2*t3^5 - 2*TH1*t3^5 - 2*TH1*t2^2*t3^3 - TH1*t2^2*t4^3 - TH1*t3^2*t4^3 - 2*TH3*t2^3*t3^2 + TH2*t3^2*t4^3 + TH3*t2^2*t4^3 + 2*TH4*t2^2*t3^3 + 2*TH4*t2^3*t3^2 + 4*TH1*t2*t3^4 + 3*TH1*t3^4*t4 - 3*TH2*t3^4*t4 - 4*TH4*t2*t3^4 + 2*TH1*t2*t3*t4^3 - 6*TH1*t2*t3^3*t4 - 2*TH3*t2*t3*t4^3 + 6*TH4*t2*t3^3*t4 + 3*TH1*t2^2*t3^2*t4 + 3*TH3*t2^2*t3^2*t4 - 6*TH4*t2^2*t3^2*t4))/(t2*(t2 - t3)*(t3 - t4)^3*(3*t2*t3 + t2*t4 - 4*t3*t4));
    
    c(3,:) = -(3*(TH2*t3^5 - TH1*t3^5 - TH1*t2^2*t3^3 - 2*TH1*t2^2*t4^3 - 2*TH1*t3^2*t4^3 + 3*TH1*t3^3*t4^2 - TH3*t2^3*t3^2 + 2*TH2*t3^2*t4^3 - 3*TH2*t3^3*t4^2 + 2*TH3*t2^2*t4^3 - TH3*t2^3*t4^2 + TH4*t2^2*t3^3 + TH4*t2^3*t3^2 + TH4*t2^3*t4^2 + 2*TH1*t2*t3^4 - 2*TH4*t2*t3^4 + 4*TH1*t2*t3*t4^3 - 4*TH3*t2*t3*t4^3 - 6*TH1*t2*t3^2*t4^2 + 3*TH1*t2^2*t3*t4^2 + 2*TH3*t2*t3^2*t4^2 + 2*TH3*t2^2*t3*t4^2 + 4*TH4*t2*t3^2*t4^2 - 5*TH4*t2^2*t3*t4^2))/(t2*(t2 - t3)*(t3 - t4)^3*(3*t2*t3 + t2*t4 - 4*t3*t4));
 
    d(3,:) = (3*TH2*t3^4 - 3*TH1*t3^4 - 3*TH1*t2^2*t3^2 - 3*TH1*t2^2*t4^2 - 3*TH1*t3^2*t4^2 - 3*TH3*t2^2*t3^2 + 3*TH2*t3^2*t4^2 + 3*TH3*t2^2*t4^2 + 6*TH4*t2^2*t3^2 + 6*TH1*t2*t3^3 + 6*TH1*t3^3*t4 - 6*TH2*t3^3*t4 - 2*TH3*t2^3*t4 - 6*TH4*t2*t3^3 + 2*TH4*t2^3*t4 + 6*TH1*t2*t3*t4^2 - 12*TH1*t2*t3^2*t4 + 6*TH1*t2^2*t3*t4 - 6*TH3*t2*t3*t4^2 + 4*TH3*t2*t3^2*t4 + 4*TH3*t2^2*t3*t4 + 8*TH4*t2*t3^2*t4 - 10*TH4*t2^2*t3*t4)/(t2*(t2 - t3)*(t3 - t4)^3*(3*t2*t3 + t2*t4 - 4*t3*t4));
 
end



