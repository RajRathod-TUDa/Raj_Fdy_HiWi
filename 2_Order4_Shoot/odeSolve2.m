function [resi, YSol] = odeSolve2(XV, bc1, bc2, bc3, bc4, lamb, outFlag)
    Y = [bc1;bc2;bc3;bc4];
    [~,b] = ode45(@odedef, XV, Y);
    resi = b(end,1)^2 + b(end,2)^2;
    if outFlag == 2
        YSol = b;
    else
        YSol = 0;
    end
    function dYdx = odedef(x,Y)
        y = Y(1);
        dydx = Y(2);
        d2ydx2 = Y(3);
        d3ydx3 = Y(4);
        dYdx(1,1) = dydx;
        dYdx(2,1) = d2ydx2;
        dYdx(3,1) = d3ydx3;
        dYdx(4,1) = (lamb/(x^4))*y;
    end
end
