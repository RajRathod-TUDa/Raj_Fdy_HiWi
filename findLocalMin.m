% This function finds local minima and returns their locations.
% Input: 2 arrays
%       X => X variable row vector (i.e. Eigenvalues)
%       Y => Y variable row vector (i.e. Residual Values)
% Output:
%       eigval => vector containing values of points in X where Y is low
% -------------------------------------------------------------------------

function eigval = findLocalMin(X,Y)
    eigval = NaN(size(X));      % Variable creation
    for i = 2:(length(Y)-1)
        if (Y(i-1)>Y(i))&&(Y(i)<Y(i+1))
            eigval(i) = X(i);       % Inserting only relevant values
        end 
    end
    eigval = rmmissing(eigval);     % filtering required values
end