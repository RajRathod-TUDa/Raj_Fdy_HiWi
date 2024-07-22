clear;
close all;
clc;

% Code for Residual values (measured at other end), for EVs taken from Collocation code.

%% Independent solution
xValues = 0:0.005:1;
ISol_1 = [0;0;1;0];
ISol_2 = [0;0;0;1];
phi = getMinors([ISol_1 ISol_2]);

%% Real and imaginary number discretization

load('eeOS.mat')        % Contains first 13 EVs from Collocation code
eeOS = eeOS(1:13,:);
cr = eeOS(:,1)';        % Real part of those EVs
ci = eeOS(:,2)';        % Imag part of those EVs

%% Solution core
% Calculates solution for all EVs and stores end residual in "resi_all"
Sol_col_1 = zeros(length(cr), length(xValues));
resi_all = zeros(1,length(cr));
for m = 1:length(cr)
        c = cr(m) + 1i*ci(m);
        [resi_all(m),Sol_col_1(m,:)] = odeSolve(xValues, phi, c,2);
end

%% Plotting
plot(cr, ci,'*b');
xlabel('real');
ylabel('imag');
grid on
text(cr,ci,'  ' + string(resi_all));
title('Residual values measured at other end, for cr & ci taken from Collocation code');

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
