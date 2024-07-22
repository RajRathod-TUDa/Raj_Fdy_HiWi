clear;
clc;

%% Independent solution
%xValues = [0 1];
xValues = 0:0.005:1;
ISol_1 = [0;0;1;0];
ISol_2 = [0;0;0;1];
phi = getMinors([ISol_1 ISol_2]);

%% Real and imaginary number taken from Colloc code
load('eeOS.mat')
eeOS = eeOS(1:13,:);
cr = eeOS(:,1)';
ci = eeOS(:,2)';

%% Solution core with Gradient optimization
Sol_col_1 = zeros(length(cr), length(xValues));
resi_all = zeros(1,length(cr));
for m = 1:1            %Currently solving for 1 EV only
        c = cr(m) + 1i*ci(m);
        [C_resi,~] = odeSolve(xValues, phi, c,1);
        LoopResiLim = 1e-16;
        LoopResi = 1;
        maxIter = 100;
        currIter = 0;
        dx = 0.02;
        dy = 0.02;
        crStore = zeros(1,maxIter);
        ciStore = zeros(1,maxIter);
    while (abs(LoopResi)>LoopResiLim) && (currIter<maxIter)
        [CRplus,~] = odeSolve(xValues, phi, c + (dx + 0*1i),2);
        [CRminus,~] = odeSolve(xValues, phi, c - (dx + 0*1i),2);
        [CIplus,~] = odeSolve(xValues, phi, c + (0 + dy*1i),2);
        [CIminus,~] = odeSolve(xValues, phi, c - (0 + dy*1i),2);

        % difference of residuals in nearby region
        dCRplus = CRplus - C_resi;
        dCRminus = CRminus - C_resi;
        dCIplus = CIplus - C_resi;
        dCIminus = CIminus - C_resi;
        dArray =  [dCRplus;dCRminus;dCIplus;dCIminus];
        % dArray is the array of residuals, in a particular order.
        % min function returns the value and tells which one is it.
        % it is used to change the value of c according to min.
        [LoopResi,chosenDiff] = min(dArray);

        if chosenDiff==1
            c = c + (dx + 0*1i);
        else
            if chosenDiff==2
                c = c - (dx + 0*1i);
            else
                if chosenDiff==3
                    c = c + (0 + dy*1i);
                else
                    c = c - (0 + dy*1i);
                end
            end
        end
        currIter = currIter + 1;

        %crstore variable stores all real values of c (used later for dplot)
        crStore(currIter) = real(c);
        ciStore(currIter) = imag(c);
    end
        [resi_all(m),Sol_col_1(m,:)] = odeSolve(xValues, phi, c,2);
end



%% Functions
function [resi, YSol]  = odeSolve(XV, Y, cComp, outFlag)
    Re = 250;
    Alp = 1;
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [~,b] = ode45(@odedef, XV, Y, options);
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
