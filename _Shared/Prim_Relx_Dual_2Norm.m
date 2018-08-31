function optmRslts = Prim_Relx_Dual_2Norm(inptStrct,optStrct,loclStrct,...
    intmRslts)
% function optmRslts = Prim_Relx_Dual(inptStrct,optStrct,loclStrct,...
%                           intmRslts)
%   (c) Yang Wang, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2018
%
% Revision: 1.0
%
% This function adopts the primal-relaxed dual global optimzation algorithm
% to solve the optimization problem of modal property difference and modal
% dynamic residual approach
%
% Input:
%   inptStrct
%     y0 - Initial values of Y variables
%     polynomials - a symbolic vector containing the constratins of the
%       model updating problem
%     numX - number of X variables
%   optStrct
%     tolGap - tolerance of the gap between lower and upper bound of the
%       objective function value
%     iterLimt - limitation on the number of iterations
%     x_lb, x_ub - lower and upper bounds of X variables
%     y_lb, y_ub - lower and upper bounds of Y variables
%     lcaSearch - whether to do local optimization
%       0 - NOT do local search
%       1 - Do local search
%     linTolbx - linear optimization toolbox
%       1 - cplex optimization toolbox from IBM
%       2 - linprog optimization toolbox from MATLAB
%     xOption - options for X variables
%       1 - use alpha as X
%       2 - use psi as X
%   loclStrct
%     lcaObj - objective function of local optimization problem
%     lcaCons - constraints of local optimization problem
%  intmRslts - structure contains intermediate results from previous runs
%     x - optimal value of X variables
%     y - optimal value of Y variables
%     nodeNum - node number of optimal result in branch and bound tree
%     numNode - total number of node in the branch and bound tree
%     iter - iteration number of primal-relaxed dual optimization algorithm
%     xRec - record of X variable values during optimization process
%     yRec - record of Y variable values during optimization process
%     treeNode - structure represents the nodes in the branch and bound tree
%     lbFval - lower bound of objective function value
%     ubFval - upper bound of objective function value
%     intRes - record intermediate results
% Output:
%   final_result - structure contains the final optimal results
%     x - optimal value of X variables
%     y - optimal value of Y variables
%     nodeNum - node number of optimal result in branch and bound tree
%     numNode - total number of node in the branch and bound tree
%     iter - iteration number of primal-relaxed dual optimization algorithm
%     xRec - record of X variable values during optimization process
%     yRec - record of Y variable values during optimization process
%     treeNode - structure represents the nodes in the branch and bound tree
%     lbFval - lower bound of objective function value
%     ubFval - upper bound of objective function value
%     intRes - store intermediate results

%% Get option values from Inpt_Strct and initialize variables for the algorithm
xIter = inptStrct.x0;
numX = length(inptStrct.x0);
numY = length(inptStrct.y0) + 1; % add mu_b to y variables

biCvx = inptStrct.polynomials(1); % biconvex polynomials - objective function
bnd = inptStrct.polynomials(2:end); % bound polynomials - X variable only

x_lb = optStrct.x_lb; x_ub = optStrct.x_ub;
y_lb = optStrct.y_lb; y_ub = optStrct.y_ub;

xOption = optStrct.xOption;
lcaSearch = optStrct.lcaSearch ;

% nonlinear solver options
nlinOptions = optStrct.nlinOptions;

lcaObj = loclStrct.lcaObj;

% Screen printer headers and figure plot
header = sprintf('\n Iter   Upp. bnd.     Low. bnd.       Iter upp. bnd.    Abs Gap    Node   Pnode\n' );
formatstr = '%5.0f   %5.3g  %13.6g    %13.6g   %13.6g  %5.0f  %5.0f\n';

h = figure;
set(h, 'WindowStyle', 'Docked');

% Symbolic for X and Y vairables
X = sym('x', [numX, 1]); Y = sym('y', [numY, 1]);

if(nargin < 4)
    % Starting from new inital value
    disp('Starting from new initial value');
    % Store intermediate results
    intRes = [];
    % Initilize lower and upper bounds of objective function value
    lbFval = 1e5;
    ubFval = inf;
    % Data structure for primal relaxed branch and bound tree
    %  -nodeLvl: level of node in the branch and bound tree
    %  -pNode: parent node in the branch and bound tree
    %  -yVal: fixed value of Y variables
    %  -linLag: linearized lagrangian function
    %  -gradLagCons: gradient of the lagrangian function
    
    % root node in the branch and bound tree
    treeNode(1) = struct('nodeLvl', 1, 'pNode', 0, 'yVal', inptStrct.y0);
    nodeNum = 1;
    numNode = 1;
    iter = 1;
    xRec = []; yRec = [];
    
else
    disp('Picking up from previous results');
    intRes = intmRslts.intRes;
    lbFval = intmRslts.lbFval;
    ubFval = intmRslts.ubFval;
    
    treeNode = intmRslts.treeNode;
    nodeNum = intmRslts.nodeNum;
    numNode = intmRslts.numNode;
    
    iter = intmRslts.iter;
    
    xRec = intmRslts.xRec; yRec = intmRslts.yRec;
    
    fvalCand = intmRslts.fval_can;
    yValCand = intmRslts.y_can;
    nodeNumCand = intmRslts.can_NodeNum;
    
end

ubLocal = inf;
fprintf(header);
ubIter = inf;

while( abs( lbFval - ubIter) > optStrct.tolGap && ...
        iter < optStrct.iterLimt)
    % NodeNum - node number of smallest lower bound value -> parent node
    % for the new nodes
    % Level for new node in the branch and bound tree
    % -> level of parent node + 1
    nLvl = treeNode(nodeNum).nodeLvl + 1;
    % Level of current node (first level to backtrack for qualifying constraints)
    tLevel = treeNode(nodeNum).nodeLvl;
    % Fixed y value for solving primal problem
    yIter = treeNode(nodeNum).yVal;
    % Construct primal problem with new fixed y value
    primProb = matlabFunction(subs(biCvx, Y(1 : end - 1), yIter),'Vars',{X});
    [xIter, primFval, ~,~, lambda] = fmincon(primProb, xIter,...
        [], [], [], [], x_lb , x_ub, [], nlinOptions);
    % Record the x and y values during optimization process
    xRec = [xRec xIter];
    yRec = [yRec [yIter; lbFval]];
    % Update upper bound of the objective function
    ubFval = min([ubFval, primFval]);
    % Do local search using updated x and y values as initial values
    if(lcaSearch)
        if(xOption == 1)
            x0 = [xIter; yIter;];
            lb = [x_lb; y_lb(1:end-1); ];
            ub = [x_ub; y_ub(1:end-1); ];
        else
            x0 = [yIter; xIter];
            lb = [y_lb(1:end-1); x_lb];
            ub = [y_ub(1:end-1); x_ub];
        end
        [localTemp, fvalTemp, localFlag] = fmincon(lcaObj, x0,...
            [], [], [], [], lb, ub, [], nlinOptions);
        % fmincon finds better local optimal value
        if(localFlag > 0 && ubLocal > fvalTemp)
            ubLocal = fvalTemp;
            xLca = localTemp;
        end
        
        % Gap between lcoal minimum and lower bound is smallet than
        % toletance value; exit the algtorithm with x, y values
        if(abs(ubLocal - lbFval) < optStrct.tolGap)
            if(xOption == 1)
                optmRslts.x = xLca(1:numX);
                optmRslts.y = xLca(numX + 1 : end-1);
            else
                optmRslts.x = xLca(numX + 1 : end-1);
                optmRslts.y = xLca(1:numX);
            end
            optmRslts.nodeNum = nodeNum;
            optmRslts.numNode = numNode;
            optmRslts.iter =  iter;
            optmRslts.xRec = xRec;
            optmRslts.yRec = yRec;
            optmRslts.treeNode = treeNode;
            optmRslts.lbFval = lbFval;
            optmRslts.ubFval = ubLocal;
            optmRslts.intRes = intRes;
            optmRslts.fvalCand = fvalCand;
            optmRslts.yValCand = yValCand;
            optmRslts.nodeNumCand = nodeNumCand;
            return
        end
        
    end
        
    % Upper bound of objective function value for current iteration
    ubIter = primFval;
    
    % Lagrange function of the new nodes
    lagFunc = biCvx + [lambda.lower;lambda.upper;]' * bnd;
    
    % Find connected variables - only X whose coeffecient contains y in the
    % Lagrangian function is connected variables
    conctX = []; xIdx = [];
    for i = 1 : numX
        % get coeffecient of x first
        [cx,tx] = coeffs(lagFunc,X(i));
        if(~isempty(tx))
            if(tx(end) == 1)
                loopX = length(cx) - 1;
            else
                loopX = length(cx);
            end
            % check whether coeffecient of x contains y 
            for j = 1:loopX
                [~,tx] = coeffs(cx(j),Y);
                if(~(length(tx) == 1 && tx(end) == 1) && ~isempty(tx))
                    conctX = [conctX X(i)];
                    xIdx = [xIdx i];
                    break;
                end
            end
        end
    end
    nonCnctIdx = setdiff(1 : numX, xIdx);
    numConctX = length(conctX);
    %% Find qualifying constraints for relaxed subproblem
    rexDualQual = sym([]);
    % Current node is the first parent node of the new node
    pNode = nodeNum;
    while(tLevel > 1)
        rexDualQual = [rexDualQual;treeNode(pNode).rexDualIter];
        tLevel = tLevel - 1;
        pNode = treeNode(pNode).pNode;
    end
    %% Binary tree for relaxed dual problem
    if(numConctX > 0)
        % number of binary tree level - numConctX
        rexDulTree = tree('root');
        for i = 1 : numConctX
            for j = 1 : (2^i)
                if(i == 1)
                    [rexDulTree, l1(j)] = rexDulTree.addnode(1, 0);
                else
                    eval(['[rexDulTree, l' num2str(i) '(j)] = rexDulTree.addnode(l' num2str(i-1) '( (ceil(j/2)) ),0);']);
                end
            end
        end
        eval(['botmNode = l' num2str(numConctX) ';']);
        
        for i = 1 : length(botmNode)
            % Construct relaxed dual problem from current iteration
            numNode = numNode + 1;
            treeNode(numNode).nodeLvl = nLvl;
            treeNode(numNode).pNode = nodeNum;
            
            treePath = findpath(rexDulTree,1,botmNode(i));
            qualX = zeros(numConctX,1);
            gradLagCons = sym(zeros(numConctX,1));
            gradLagFunc = sym(zeros(1,numConctX));
            % find value of Conncected X variables
            % find relaxed Lagrangian function and companying constraints
            for j = 1 : numConctX
                % differentiate Lagrangian function w.r.t connected variable
                gradLagFunc(1,j) = subs(diff(lagFunc,conctX(j)), X, xIter);
                if( rem(treePath(1 + j ),2) == 0 )
                    gradLagCons(j,1) = gradLagFunc(1,j);
                    qualX(j) = x_ub(xIdx(j));
                else
                    gradLagCons(j,1) = -gradLagFunc(1,j);
                    qualX(j) = x_lb(xIdx(j));
                end
                
            end
            %    lagFunc @ x_K + graident of lagFunc * (x_b - xk) - mu_b
            linLag = subs(lagFunc,X,xIter) + gradLagFunc * (qualX - xIter(xIdx)) - Y(end) ;
            
            rexDualIter = [linLag; gradLagCons];
            
            treeNode(numNode).rexDualIter = rexDualIter;
                        
            if(nLvl == 2)
                rexDual = rexDualIter;
            else
                rexDual = [rexDualQual; rexDualIter];
            end
            ceq = 0;
            
            accpCons = matlabFunction(rexDual,ceq,'Outputs',{'cineq','ceq'},'Vars',{Y});
            % minimize mu_B
            objDual = matlabFunction(Y(end), 'Vars', {Y});
            
            [yTemp, ~, exitflag] = fmincon(objDual,...
                [yIter; 0], [], [], [], [], y_lb, y_ub, accpCons, nlinOptions);
                        
            if(exitflag > 0)
                % save y and lower bound value
                treeNode(numNode).yVal = yTemp(1 : numY - 1);
                treeNode(numNode).fval = yTemp(end);
                % ToDO: how useful is following code
                % make sure y value is updated from the iteration
                if(norm(yTemp(1 : numY - 1) - yIter) > 1e-7)
                    intRes = [intRes [yTemp; numNode]];
                end
            end
        end
    else
        % No conncected X variables
        numNode = numNode + 1;
        treeNode(numNode).nodeLvl = nLvl;
        treeNode(numNode).pNode = nodeNum;
        
        linLag = subs(lagFunc,X,xIter) - Y(end) ;
        rexDualIter = linLag;
        
        treeNode(numNode).rexDualIter = rexDualIter;
                
        if(nLvl == 2)
            rexDual = rexDualIter;
        else
            rexDual = [rexDualQual;rexDualIter];
        end
        
        accpCons = matlabFunction(rexDual,ceq,'Outputs',{'cineq','ceq'},'Vars',{Y});
        
        % minimize mu_B
        objDual = matlabFunction(Y(end),'Vars',{Y});
        
        [yTemp, ~, exitflag] = fmincon(objDual,...
            [yIter;0], [], [], [], [],y_lb , y_ub, accpCons, nlinOptions);
        
        
        if(exitflag > 0)
            % save y and lower bound value
            treeNode(numNode).yVal = yTemp(1 : numY - 1);
            treeNode(numNode).fval = yTemp(end);
            % ToDO: how useful is following code
            % make sure y value is updated from the iteration
            if(norm(yTemp(1 : numY - 1) - yIter) > 1e-7)
                intRes = [intRes [yTemp; numNode]];
            end
        end
        
    end
    
    % Fathom node with lb value larger than ub
    fathIdx = intRes(numY,:) > (min(ubFval,ubLocal) + optStrct.tolGap);
    intRes(:,fathIdx) = [];
    
    % Find new lower bound value
    lbFvalIter = min( intRes( numY, :));
    
    % Find candidate with same lower bound value
    lbFvalIdx = find(abs(intRes(numY,:) - lbFvalIter) < 1e-7);
    
    % Pick y candidates from current iteration achieve the same lower
    % bound value
    yValIter = intRes(1 : numY - 1, lbFvalIdx);
    % Pick NodeNum of y candidates from current iteration
    nodeNumIter = intRes(end, lbFvalIdx);
    % Substract old candidate from last iteration if same lower bound
    % is achieved at current iteration
    if(abs(lbFvalIter - lbFval) < 1e-7)
        % Find duplicate index bewtween new and old y candidate
        [~,redIdx,~] = intersect(nodeNumIter', nodeNumCand', 'rows');
        remIdx = setdiff(1:length(nodeNumIter),redIdx,'stable');
        nodeNumIter = nodeNumIter(remIdx);
        yValIter = yValIter(:,remIdx);
    else
        % Construct new candidate set if new lower bound value is achieved
        % at current iteration
        fvalCand = [];
        yValCand = [];
        nodeNumCand = [];
        lbFval = lbFvalIter;
    end
    
    % Calculate upper bound from all new y candidates
    fvalCanIter = zeros(1,size(yValIter,2));
    
    for i = 1:size(yValIter,2)
        primProb = matlabFunction(subs(biCvx, Y(1 : end - 1), yValIter(:,i)),'Vars',{X});
        [~, fvalCanIter(i)] = fmincon(primProb, xIter,...
            [], [], [], [], x_lb , x_ub,[],nlinOptions);
    end
    % Add upper bound value from new y candidate to candidate group
    fvalCand = [fvalCand fvalCanIter];
    yValCand = [yValCand yValIter];
    nodeNumCand = [nodeNumCand nodeNumIter];
    
    % Find smallest upper bound from all y (including current and previous
    % iterations) to be used for next iteration
    [~,min_idx] = min(fvalCand);
    nodeNum = nodeNumCand(min_idx);
    
    % Remove next y value from candidate group
    fvalCand(min_idx) = [];
    yValCand(:,min_idx) = [];
    nodeNumCand(min_idx) = [];
    
    % Remove next y value from all result group
    intRes(:,intRes(end,:) == nodeNum) = [];
    
    fprintf(formatstr,iter,min(ubFval,ubLocal),lbFval,ubIter,ubIter - lbFval,nodeNum,treeNode(nodeNum).pNode );
    iter = iter + 1;
    stem(1:numX,xIter(1:numX),'o')
    ylim([min(x_lb(1:numX));max(x_ub(1:numX))]);
    xticks(1:1:numX);
    pause(0.5);
end
optmRslts.x = xIter;
optmRslts.y = yIter;
optmRslts.nodeNum = nodeNum;
optmRslts.numNode = numNode;
optmRslts.iter =  iter;
optmRslts.xRec = xRec;
optmRslts.yRec = yRec;
optmRslts.treeNode = treeNode;
optmRslts.lbFval = lbFval;
optmRslts.ubFval = ubFval;
optmRslts.intRes = intRes;
optmRslts.fvalCand = fvalCand;
optmRslts.yValCand = yValCand;
optmRslts.nodeNumCand = nodeNumCand;
end

