function [final_L] = OuterSecMeth(XV, bc1, bc2, bc3, bc4, L_St, L_End, bc5, bc6)

    L_curr = L_St;
    tempFunc = @(abc) odeSolve(XV, bc1, bc2, abc, bc4, L_curr, 1);
    yDD_St_corr = SimpleSecMeth(tempFunc, bc3, 0);
    [~,fullSol] = odeSolve(XV, bc1, bc2, yDD_St_corr, bc4, L_curr, 2);
    resi1 = abs(fullSol(end,1) - bc5);
    resi2 = abs(fullSol(end,2) - bc6);
    F1 = resi1 + resi2;

    L_curr = L_End;
    tempFunc = @(abc) odeSolve(XV, bc1, bc2, abc, bc4, L_curr, 1);
    yDD_St_corr = SimpleSecMeth(tempFunc, bc3, 0);
    [~,fullSol] = odeSolve(XV, bc1, bc2, yDD_St_corr, bc4, L_curr, 2);
    resi1 = abs(fullSol(end,1) - bc5);
    resi2 = abs(fullSol(end,2) - bc6);
    F2 = resi1 + resi2;

    funcResi = 1e-12;
    maxIter = 200;
    currIter = 0;

    L1 = L_St;
    L2 = L_End;

    while ((abs(F2-F1))>funcResi) && (currIter<maxIter)
        L3 = L2 - (F2)*((L2 - L1)/(F2 - F1));
        F1 = F2;

        L_curr = L3;
        tempFunc = @(abc) odeSolve(XV, bc1, bc2, abc, bc4, L_curr, 1);
        yDD_St_corr = SimpleSecMeth(tempFunc, bc3, 0);
        [~,fullSol] = odeSolve(XV, bc1, bc2, yDD_St_corr, bc4, L_curr, 2);
        resi1 = abs(fullSol(end,1) - bc5);
        resi2 = abs(fullSol(end,2) - bc6);
        F2 = resi1 + resi2;

        L1 = L2;
        L2 = L3;
        currIter = currIter + 1
    end
    final_L = L2;
end