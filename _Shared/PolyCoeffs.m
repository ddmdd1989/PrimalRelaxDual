function [coeffX, constX, coeffY, constY] = PolyCoeffs(polynomial, numX, numY )

% function [coeffX, constX, coeffY, constY] = PolyCoeffs(polynomial, numX, numY )
%
%   (c) Yang Wang, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2018
%
% Revision: 1.0
%
% This function obtains the coeffecients of the X and Y variables, and constants
% as well as of the polynomial
%
% Input:
%   polynomial - a symbolic vector (numPoly x 1) contains two groups of
%   variables, X and Y;
%   numX - number of X variables
%   numY - number of Y variables
%
% Output:
%   coeffX - coeffecient of X variables (numPoly x (numY + 1) x numX)
%      Dimension 1: corresponding to each polynomial
%      Dimension 2: 1 -> numY - coeffecient of biliear terms
%                   numY + 1 - coeffecient of X linear term
%      Dimension 3: corresponding to each X variable
%   constX - terms not containing X variables (numPoly x (numY + 1))
%       Dimension 1: corresponding to each polynomial
%       Dimension 2: 1 -> numY - coeffecient of Y linear terms
%                    numY + 1 - constant term
%   coeffY - coeffecient of Y variables (numPoly x numX x numY)
%      Dimension 1: corresponding to each polynomial
%      Dimension 2: 1 -> numX - coeffecient of biliear terms 
%                   numX + 1 - coeffecient of Y linear term
%      Dimension 3: corresponding to each Y variable
%   constY - terms not containing Y variables (numPoly x (numX + 1))
%       Dimension 1: corresponding to each polynomial
%       Dimension 2: 1 -> numX - coeffecient of X linear terms
%                    numX + 1 - constant term
X = sym('x',[numX, 1]);
Y = sym('y',[numY, 1]);

numPoly = length(polynomial);

coeffX = zeros(numPoly, numY + 1, numX);
constX = zeros(numPoly, numY + 1);

for j = 1 : numPoly
    % Find coeffecient for each X variable
    [cx,tx] = coeffs(polynomial(j),X);
    % if last entry in tx is 1 -> polynomial has terms not containing X
    if(tx(end) == 1)
        % Find coeffect of Y variable
        [cy,ty] = coeffs(cx(end), Y);
        [~,ia,ib] = intersect(ty,Y);
        constX(j, ib) = double(cy(ia));
        if(ty(end) == 1)
            constX(j, end) = double(cy(end));
        end
        [~, ~, biLinIdx] = intersect(tx(1 : end - 1), X);
    else
        [~, ~, biLinIdx] = intersect(tx, X);
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
% 
% coeffY = zeros(numPoly, numX, numY);
% constY = zeros(numPoly, numX, 1);


coeffY = zeros(numPoly, numX + 1, numY);
constY = zeros(numPoly, numX + 1, 1);

for j = 1 : numPoly
    [cy,ty] = coeffs(polynomial(j), Y);
    if(ty(end) == 1)
        [cx,tx] = coeffs(cy(end), X);
        [~,ia,ib] = intersect(tx, X);
        constY(j, ib) = double(cx(ia));
        if(tx(end) == 1)
            constY(j, end) = double(cx(end));
        end
        [~, ~, biLinIdx] = intersect(ty(1 : end - 1), Y);
    else
        [~, ~, biLinIdx] = intersect(ty(1 : end), Y);
    end
    for i = 1 : length(biLinIdx)
        [cx,tx] = coeffs(cy(i), X);
        [~,ia,ib] = intersect(tx, X);
        coeffY(j, ib, biLinIdx(i)) = double(cx(ia));
        if(tx(end) == 1)
            coeffY(j, end, biLinIdx(i)) = double(cx(end));
        end
    end
end
end

