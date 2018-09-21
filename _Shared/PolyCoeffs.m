function [coeffX, constX] = PolyCoeffs(polynomial, numX, numY)

% function [coeffX, constX, coeffY, constY] = PolyCoeffs(polynomial, numX, numY )
%
%   (c) Yang Wang, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2018
%
% Revision: 1.0
%
% For a bilinear function whose bilinear terms are only between an X_i and
% a Y_j,  this MATLAB routine extracts the coeffecients of the X and Y
% variables, and constants of the polynomials
%
% Input:
%   polynomial - a symbolic vector (numPoly x 1) contains two groups of
%   variables, X and Y;
%   numX - number of X variables
%   numY - number of Y variables
%
% Output:
%   coeffX - coeffecient of X variables (numPoly x (numY + 1) x numX)
%      1st dimension (numPoly): corresponding to all polynomials
%      3rd dimension (numX): corresponding to all X variables in each polynomial
%      2nd dimension (numY + 1): for each X_i variable in each polynomial,
%           1 ~ numY: coeffecients of bilinear terms involving X_i * Y_j
%           numY + 1: coeffecient of the linear term involing only X_i
%   constX - terms not containing X variables (numPoly x (numY + 1))
%       1st dimension (numPoly): corresponding to all polynomials
%       2nd dimension (numY + 1): for each polynomial,
%           1 ~ numY: coeffecient of linear terms involving only Y_j 
%           numY + 1: constant term that doesn't contain any X_i or Y_j
 
X = sym('x',[numX, 1]);
Y = sym('y',[numY, 1]);

numPoly = length(polynomial);

coeffX = zeros(numPoly, numY + 1, numX);
constX = zeros(numPoly, numY + 1);

for j = 1 : numPoly
    % Find coeffecient for each X variable
    [cx,tx] = coeffs(polynomial(j),X);
    % if last entry in tx is 1 -> polynomial has terms not containing X
    % -> only linear in y and constant only
    if(tx(end) == 1)
        % Find coeffect of Y variable
        [cy,ty] = coeffs(cx(end), Y);
        [~,ia,ib] = intersect(ty,Y);
        constX(j, ib) = double(cy(ia));
        if(ty(end) == 1)
            constX(j, end) = double(cy(end));
        end
        [~, ~, biLinIdx] = intersect(tx(1 : end - 1), X, 'stable');
    else
        [~, ~, biLinIdx] = intersect(tx, X, 'stable');
    end
    for i = 1 : length(biLinIdx)
        % coeffecient for each Y variable
        [cy,ty] = coeffs(cx(i), Y);
        [~,ia,ib] = intersect(ty,Y);
        coeffX(j, ib, biLinIdx(i)) = double(cy(ia));
        if(ty(end) == 1)
            coeffX(j, end, biLinIdx(i)) = double(cy(end));
        end
    end
end
end

