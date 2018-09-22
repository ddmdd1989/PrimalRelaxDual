if(nLvl == 2)
    % At second level, the root node is y0 without qualifying
    % constraints
    ADual = ADualK;
    bDual = bDualK;
else
    ADual = [AQual; ADualK];
    bDual = [bQual; bDualK];
end

% Solve Problem 2:  minimize c' * {y; muB} = muB
c = [zeros(numY,1);
    1];
if(optimzOpts.linTolbx == 1)
    [yTemp, ~, exitflag] = cplexlp(c, ADual, bDual, [], [],...
        y_lb, y_ub, [], linOptions);
else
    [yTemp, ~, exitflag] = linprog(c, ADual, bDual, [], [],...
        y_lb, y_ub, [], linOptions);
end
