% Fixed yK value for solving primal problem
yK = treeNode(nodeNum).yVal;
% Construct primal problem with new fixed yK value
%
%       minimize  c' * x
%          x
%     subject to  A * x <= b
%            x_lb <= x <= x_ub
A_prim = zeros(size(coeffX,1),numX + 1);
for i = 1: numX + 1
    A_prim(:,i) = coeffX(:,:,i) * [yK; 1];
end
b_prim = -constX * [yK; 1];
c = [zeros(numX,1); 1];
% Solve primal optimization problem
% lambda:
% lambda.ineqlin: muK
if(optimzOpts.linTolbx == 1)
    [xPrim, ubFvalK, ~, ~, lambda] = cplexlp(c, A_prim, b_prim, [], [],...
        x_lb, x_ub, [], linOptions);
else
    [xPrim, ubFvalK, ~, ~, lambda] = linprog(c, A_prim, b_prim, [], [],...
        x_lb, x_ub, [], linOptions);
end