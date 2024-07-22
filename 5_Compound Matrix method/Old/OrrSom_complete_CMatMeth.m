clear;
close all;
clc;

xValues = [0 1];
ISol_1 = [0;0;1;0];
ISol_2 = [0;0;0;1];
phi = getMinors([ISol_1 ISol_2]);

ci_St = 0.003;
ci_end = 0.004;

tempFunc = @(abc) odeSolve(xValues, phi, abc);

ci = ci_St: 0.0001 : ci_end;
for l = 1:length(ci)
    s(l) = odeSolve(xValues, phi, 0.2375263 + 1i*ci(l));
end
plot(ci,log(abs(s)));

function resi = odeSolve(XV, Y, cComp)
    Re = 10000;
    Alp = 1;
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [solx,b] = ode45(@odedef, XV, Y, options);
    resi = abs(b(end,1))^2;
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
