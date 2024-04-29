function [final_x] = SimpleSecMeth(funcHand,initGs1,initGs2)

    funcResi = 1e-16;
    F1 = funcHand(initGs1);
    F2 = funcHand(initGs2);
    maxIter = 30;
    currIter = 1;
    x1 = initGs1;
    x2 = initGs2;

    while ((abs(F2-F1))>funcResi) && (currIter<maxIter)
        x3 = x2 - (F2)*((x2 - x1)/(F2 - F1));
        F1 = F2;
        F2 = funcHand(x3);
        x1 = x2;
        x2 = x3;
        currIter = currIter + 1;
    end
    final_x = x2;
end