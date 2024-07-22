clear;
close all;
clc;

%% Independent solution
%xValues = [0 1];
xValues = 0:0.005:1;
ISol_1 = [0;0;1;0];
ISol_2 = [0;0;0;1];
phi = getMinors([ISol_1 ISol_2]);

%% Real and imaginary number discretization
cr = -1:0.05:2;
ci = -8:0.05:-4;

%% Solution core
Sol_col_1 = zeros(length(cr),length(ci),length(xValues));
for m = 1:length(cr)
    for n = 1:length(ci)
        c = cr(m) + 1i*ci(n);
        [~,Sol_col_1(m,n,:)] = odeSolve(xValues, phi, c,2);
    end
end


Sol_col_1 = log(abs(Sol_col_1));
last = Sol_col_1(:,:,201);

%% Plotting
surface(ci,cr,last);
view(3);
colormap;
xlabel('ci');
ylabel('cr');
title('log(abs(residual values at other end))')

%% Functions
function [resi, YSol]  = odeSolve(XV, Y, cComp, outFlag)
    Re = 250;
    Alp = 1;
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [solx,b] = ode45(@odedef, XV, Y, options);
    resi = abs(b(end,1))^2;
    if outFlag == 2
        YSol = b(:,1);
    else
        YSol = 0;
    end
    
    function B = odedef(x,Y)
        T1 = (2*(Alp^2) + Re*1i*Alp*((2*x-x^2)-cComp));
        T2 = -(Alp^4 + Alp*Re*1i*((Alp^2)*(2*x-x^2 - cComp) - 1));
        B(1) = Y(2);
        B(2) = Y(3) + Y(4);
        B(3) = T1*Y(2) + Y(5);
        B(4) = Y(5);
        B(5) = -T2*Y(1) + T1*Y(4) + Y(6);
        B(6) = -T2*Y(2);
        B=B';
    end
end
