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
px2 = px1 - 30e-3;
py2 = py1;
pz2 = pz1 + 45e-3;

% 2ND INTERMEDIATE POINT 
px3 = px4;
py3 = py1 - 20e-3;
pz3 = pz4 + 30e-3;

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
[A,B,C,D] = tranjectory(M1,M2,M3,M4,t2,t3,t4);

n = 1e2;

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

%% TRANJECTORY 
% CALCULATION OF WORKING EDGE TRANJECTORY  
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


%% PLOT 3D TRANJECTORY 
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
title("Working Edge Tranjectory -  Bottle 1");
xlabel("x (m)");
ylabel("y (m)");
zlabel("z (m)");
grid on 

%% CREATE VIDEO WHITH ROBOT TRANJECTORY 

iloc = [20e-3 0 0 1]';
jloc = [0 20e-3 0 1]';
kloc = [0 0 20e-3 1]';

iloc1 = [30e-3 0 0 1]';
jloc1 = [0 30e-3 0 1]';
kloc1 = [0 0 30e-3 1]';

iloc2 = [25e-3 0 0 1]';
jloc2 = [0 25e-3 0 1]';
kloc2 = [0 0 25e-3 1]';

v = VideoWriter('Tranjectory_for_bottle 1.avi');
v.Quality = 100;
v.FrameRate = 60;
open(v);

figure('Name','Tranjectory for bottle 1');

for cnt = 1:n
    plot3(x,y,z,"LineWidth",1.5);
    hold on 
    Th = zeros(n,6);
    % 1ST POLYNOMIAL 
    if (T(cnt) <= t2)
        for j = 1:6
            Th(cnt,j) = A(1,j) + B(1,j)*T(cnt) + C(1,j)*T(cnt)^2 + D(1,j)*T(cnt)^3;
        end
        
    elseif (T(cnt) <= t3) && (T(cnt) > t2)
        %2ND POLYNOMIAL 
       for j = 1:6
            Th(cnt,j) = A(2,j) + B(2,j)*T(cnt) + C(2,j)*T(cnt)^2 + D(2,j)*T(cnt)^3;
        end
      
    elseif (T(cnt) <= t4) && (T(cnt) > t3)
        %2ND POLYNOMIAL 
        for j = 1:6
            Th(cnt,j) = A(3,j) + B(3,j)*T(cnt) + C(3,j)*T(cnt)^2 + D(3,j)*T(cnt)^3;
        end
    end
    
    [M01,M12,M23,M34,M45,M56] = direct_kinematic(Th(cnt,1),Th(cnt,2),Th(cnt,3),Th(cnt,4),Th(cnt,5),Th(cnt,6));
    
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
    
    plot3([0 20e-3],[0 0],[0 0],"Color","blue","LineWidth",2);
    hold on
    text(20e-3,0,0,"x_0");
    hold on
    plot3([0 0],[0 20e-3],[0 0],"Color","green","LineWidth",2);
    hold on
    text(0,20e-3,0,"y_0");
    hold on 
    plot3([0 0],[0 0],[0 20e-3],"Color","red","LineWidth",2);
    hold on 
    text(0,0,20e-3,"z_0");
    hold on 
    plot3([0 P1(1)],[0 P1(2)],[0 P1(3)]);
    hold on
    plot3([P1(1) i1(1)],[P1(2) i1(2)],[P1(3) i1(3)],"Color","blue","LineWidth",2);
    hold on
    text(i1(1),i1(2),i1(3),"x_1");
    hold on 
    plot3([P1(1) j1(1)],[P1(2) j1(2)],[P1(3) j1(3)],"Color","green","LineWidth",2);
    text(j1(1),j1(2),j1(3),"y_1");
    hold on 
    plot3([P1(1) k1(1)],[P1(2) k1(2)],[P1(3) k1(3)],"Color","red","LineWidth",2);
    hold on
    text(k1(1),k1(2),k1(3),"z_1");
    hold on
    
    plot3([P1(1) P2(1)],[P1(2) P2(2)],[P1(3) P2(3)]);
    hold on 
    plot3([P2(1) i2(1)],[P2(2) i2(2)],[P2(3) i2(3)],"Color","blue","LineWidth",2);
    hold on 
    text(i2(1),i2(2),i2(3),"x_2");
    hold on 
    plot3([P2(1) j2(1)],[P2(2) j2(2)],[P2(3) j2(3)],"Color","green","LineWidth",2);
    hold on
    text(j2(1),j2(2),j2(3),"y_2");
    hold on
    plot3([P2(1) k2(1)],[P2(2) k2(2)],[P2(3) k2(3)],"Color","red","LineWidth",2);
    hold on
    text(k2(1),k2(2),k2(3),"z_2");
    hold on

    plot3([P2(1) P3(1)],[P2(2) P3(2)],[P2(3) P3(3)]);
    hold on

    plot3([P3(1) i3(1)],[P3(2) i3(2)],[P3(3) i3(3)],"Color","blue","LineWidth",2);
    hold on
    text(i3(1),i3(2),i3(3),"x_3");
    hold on
    plot3([P3(1) j3(1)],[P3(2) j3(2)],[P3(3) j3(3)],"Color","green","LineWidth",2);
    hold on
    text(j3(1),j3(2),j3(3),"y_3");
    hold on 
    plot3([P3(1) k3(1)],[P3(2) k3(2)],[P3(3) k3(3)],"Color","red","LineWidth",2);
    hold on 
    text(k3(1),k3(2),k3(3),"z_3");
    hold on

    p3 = [30e-3 0 0 1]';
    P34 = M01*M12*M23*p3;
   

    plot3([P3(1) P34(1)],[P3(2) P34(2)],[P3(3) P34(3)]);
    hold on
    plot3([P34(1) P4(1)],[P34(2) P4(2)],[P34(3) P4(3)]);
    hold on

    plot3([P4(1) i4(1)],[P4(2) i4(2)],[P4(3) i4(3)],"Color","blue","LineWidth",2);
    hold on
    text(i4(1),i4(2),i4(3),"x_4");
    hold on
    plot3([P4(1) j4(1)],[P4(2) j4(2)],[P4(3) j4(3)],"Color","green","LineWidth",2);
    hold on
    text(j4(1),j4(2),j4(3),"y_4");
    hold on
    plot3([P4(1) k4(1)],[P4(2) k4(2)],[P4(3) k4(3)],"Color","red","LineWidth",2);
    hold on
    text(k4(1),k4(2),k4(3),"z_4");
    hold on 

    plot3([P4(1) P5(1)],[P4(2) P5(2)],[P4(3) P5(3)]);
    hold on

    plot3([P5(1) i5(1)],[P5(2) i5(2)],[P5(3) i5(3)],"Color","blue","LineWidth",2);
    hold on
    text(i5(1),i5(2),i5(3),"x_5");
    hold on
    plot3([P5(1) j5(1)],[P5(2) j5(2)],[P5(3) j5(3)],"Color","green","LineWidth",2);
    hold on 
    text(j5(1),j5(2),j5(3),"y_5");
    hold on 
    plot3([P5(1) k5(1)],[P5(2) k5(2)],[P5(3) k5(3)],"Color","red","LineWidth",2);
    hold on 
    text(k5(1),k5(2),k5(3),"z_5");
    hold on 

    plot3([P5(1) P6(1)],[P5(2) P6(2)],[P5(3) P6(3)]);
    hold on 
    
    plot3([P6(1) i6(1)],[P6(2) i6(2)],[P6(3) i6(3)],"Color","blue","LineWidth",2);
    hold on 
    text(i6(1),i6(2),i6(3),"x_6");
    hold on 
    plot3([P6(1) j6(1)],[P6(2) j6(2)],[P6(3) j6(3)],"Color","green","LineWidth",2);
    hold on 
    text(j6(1),j6(2),j6(3),"y_6");
    hold on 
    plot3([P6(1) k6(1)],[P6(2) k6(2)],[P6(3) k6(3)],"Color","red","LineWidth",2);
    hold on 
   
    text(k6(1),k6(2),k6(3),"z_6");
 
    hold on 
    
    view(3);
    
    title("Tranjectory for bottle 1");
    xlabel("x (m)");
    ylabel("y (m)");
    zlabel("z (m)");
    axis([-0.15 0.3 -0.15 0.3 0 0.45]); 
    grid on 
    
    hold off
    
    pause(0.1);
    
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);

%% FUNCTION FOR THE SOLUTION OF THE DIRECT KINEMATIC 

function [M01,M12,M23,M34,M45,M56] = direct_kinematic(th1,th2,th3,th4,th5,th6)
    % ROBOT D-H PARAMETERS
    a = ([0 -90 0 -90 90 -90])*pi/180;
    L = ([0 0 210 30 0 0])*(1e-3);
    d = ([80 + 103 0 0 180 + 41.5 0 23.7])*(1e-3);
    
    th = [th1 (th2-pi/2) th3 th4 th5 th6];
    
    M = zeros(4,4,6);
    
    % TRANSFORMS 
    for j = 1:6
         M(:,:,j) = [cos(th(j)) -sin(th(j)) 0 L(j) ; sin(th(j))*cos(a(j)) cos(th(j))*cos(a(j)) -sin(a(j)) -sin(a(j))*d(j) ; sin(th(j))*sin(a(j)) cos(th(j))*sin(a(j)) cos(a(j)) cos(a(j))*d(j) ; 0 0 0 1];
    end
    
    M01 = M(:,:,1);
    M12 = M(:,:,2);
    M23 = M(:,:,3);
    M34 = M(:,:,4);
    M45 = M(:,:,5);
    M56 = M(:,:,6);
end
%% FUNCTION FOR THE CALCULATION OF POLYNOMIAL INTERPPOLATION COEFFICIENTS 
% CALCULATION OF POLYNOMIAL COEFFICIENTS 
function [a,b,c,d] = tranjectory(M1,M2,M3,M4,t2,t3,t4)
    
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

