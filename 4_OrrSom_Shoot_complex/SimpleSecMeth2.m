function [final_x,F2_vals] = SimpleSecMeth2(funcHand,initGs1,initGs2)

    funcResi = 1e-16;
    F1 = funcHand(initGs1);
    F2 = funcHand(initGs2);
    maxIter = 200;
    currIter = 0;
    x1 = initGs1;
    x2 = initGs2;
    F2_vals = NaN(1,maxIter);

    while ((abs(F2-F1))>funcResi) && (currIter<maxIter)
        x3 = x2 - (F2)*((x2 - x1)/(F2 - F1));
        F1 = F2;
        F2 = funcHand(x3);
        x1 = x2;
        x2 = x3;
        currIter = currIter + 1;
        F2_vals(currIter) = F2;
    end
    final_x = x2;
end

