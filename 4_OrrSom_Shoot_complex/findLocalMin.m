% This function finds local minima and returns their locations.
% Input: 2 arrays
%       X => row vector (i.e. Eigenvalues)
%       Y => row vector (i.e. Logarithm of Residual Values)
% Output:
%       eigval => vector containing values of points in X, where Y is low.
%       indxEV => index where above value was found in X.
%
% Working:
% --1) create variable "eigval" and initialize with NaN
% --2) j is the number of neighbouring points the code has to consider for
%   further calculations.
% --3) Run in loop from 2nd value to second-last value. Eigenvalues will
%   fulfill following conditions:
%
%   1st IF: check if adjacent point values are bigger than current.
%       This is the location of a local minima.
%   Note: To differentiate between 'sharp dip', 'Slow dip' or oscillation
%   minimas, the inner 'if blocks' are present.
%
%   2nd IF: for this it calculates difference between residual values of
%   'j' points before and after current point. Every point must be smaller
%   than its previous in (i-j to i) region and similar logic for (i to i+j)
%   region. This confirms that it was continous downward trend before
%   minima and continous upward trend after minima. This usually restricts
%   values with local oscillations. If required, the value of j can be
%   increased.
%
%   3rd IF: only if the difference between the residual values of
%   (i-j), i and (i+j) points is more than a 1 magnitude
%   (i.e. 10^1 if no log is taken), the i th point is a sharp minima.
% -------------------------------------------------------------------------

function [eigval,indxEV] = findLocalMin(X,Y)
    eigval = NaN(size(X));      % Variable creation
    indxEV = NaN(size(X));      % Variable creation
    j = 10;
    for i = 2:(length(Y)-1)
        if (Y(i-1)>Y(i))&&(Y(i)<Y(i+1))
            diffVal = diff(Y((i-j):(i+j)));
            if (all(diffVal(1:j)<0)) && (all(diffVal(j+1:end)>0))
                if (abs(Y(i)-Y(i-j))>1 && abs(Y(i)-Y(i+j))>1)
                    eigval(i) = X(i);   % Inserting only relevant values
                    indxEV(i) = i;
                end
            end
        end 
    end
    eigval = rmmissing(eigval);
    indxEV = rmmissing(indxEV);
end