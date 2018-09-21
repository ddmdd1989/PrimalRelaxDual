function optmRslts = Prim_Relx_Dual(optimzProb, optimzOpts, localOptm,...
    intmRslts)
% function optmRslts = Prim_Relx_Dual(optimzProb,optimzOpts,localOptm,...
%                           intmRslts)
%   (c) Yang Wang, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2018
%
% Revision: 1.0
%
% This function adopts the primal-relaxed dual global optimization algorithm
% to solve the optimization problem of modal property difference and modal
% dynamic residual approach
%
% Input:
%   optimzProb
%     y0 - Initial values of Y variables
%     Information of constraints
%       coeffX - coeffecient of X variables (numPoly x (numY + 1) x numX)
%         Row: corresponding to each constaint
%         Col: 1 ~ numY - coeffecient of biliear terms
%                  numY + 1 - coeffecient of X linear term
%         Dimension 3: corresponding to each X variable
%       constX - terms not containing X variables (numPoly x (numY + 1))
%                linear in Y or constaints
%         Row: corresponding to each constraint
%         Col: 1 ~ numY - coeffecient of Y linear terms
%              numY + 1 - constant term
%   optimzOpts
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
%   localOptm
%     lcaObj - objective function of original optimization problem
%     lcaCons - constraints of original optimization problem
%   intmRslts - structure contains intermediate results from previous runs
%     x - optimal value of X variables
%     yK - optimal value of Y variables
%     nodeNum - node number of optimal result in branch and bound tree
%     numNodes - total number of node in the branch and bound tree
%     K - iteration number of primal-relaxed dual optimization algorithm
%     xRec - record of X variable values during optimization process
%     yRec - record of Y variable values during optimization process
%     treeNode - structure represents the nodes in the branch and bound tree
%     lbFval - lower bound of objective function value
%     ubFval - upper bound of objective function value
%     intRes - record intermediate results
% Output:
%   optmRslts - structure contains the final optimal result
%     x - optimal value of X variables
%     yK - optimal value of Y variables
%     nodeNum - node number of optimal result in branch and bound tree
%     numNodes - total number of node in the branch and bound tree
%     K - iteration number of primal-relaxed dual optimization algorithm
%     xRec - record of X variable values during optimization process
%     yRec - record of Y variable values during optimization process
%     treeNode - structure represents the nodes in the branch and bound tree
%     lbFval - lower bound of objective function value
%     ubFval - upper bound of objective function value
%     intRes - store intermediate results

%% Get option values from Inpt_Strct and initialize variables for the algorithm

x_lb = optimzOpts.x_lb; x_ub = optimzOpts.x_ub;
y_lb = optimzOpts.y_lb; y_ub = optimzOpts.y_ub;

numX = length(x_lb) - 1; % excluding the delta term
numY = length(optimzProb.y0);

coeffX = optimzProb.coeffX;
constX = optimzProb.constX;

% Screen printer headers and figure plot
header = sprintf('\n Iter   Upp. bnd.     Low. bnd.       Iter upp. bnd.    Abs Gap    Node   Pnode\n' );
formatstr = '%5.0f   %5.6g  %13.6g    %13.6g   %13.6g  %5.0f  %5.0f\n';

h = figure;
set(h, 'WindowStyle', 'Docked');

% linear solver option
if (optimzOpts.linTolbx == 1)
    % Cplex options
    linOptions = cplexoptimset('cplex');
    linOptions = cplexoptimset(linOptions,'lpmethod',1);
else
    % linprog options
    linOptions = optimoptions('linprog', 'Display', 'off',...
        'MaxIterations', 1e8, 'ConstraintTolerance', 1e-9);
end
% nonlinear solver options - do local search with fmincon
nlinOptions = optimoptions('fmincon','Display','off',...
    'algorithm','sqp','MaxFunEvals',50000,'MaxIter',10000);

if(nargin < 4)
    % Starting from new inital value
    disp("Starting from new initial value");
    % Store intermediate results
    intRes = [];
    % Initilize lower and upper bounds of objective function value
    lbFval = 1e5;
    ubFval = inf;
    % Data structure for primal relaxed branch and bound tree
    %  -nodeLvl: level of node in the branch and bound tree
    %  -pNode: parent node in the branch and bound tree
    %  -yVal: fixed value of Y variables
    % For the optimization problem minimizing muB, all affine inequalities
    % at the current K-th iteration, including the qualifying Lagrangian
    % and accompanying derivative constraints while excluding the bounds,
    % are expressed as
    %            [ADualK] * {y; muB} <= {bDualK}
    %  -ADualK : matrix [ADualK]
    %  -bDualK : vector {bDualK}
    
    % root node in the branch and bound tree
    treeNode(1) = struct('nodeLvl', 1, 'pNode', 0, 'yVal', optimzProb.y0,...
        'ADualK', [], 'bDualK', []);
    nodeNum = 1;
    numNodes = 1;
    K = 0;
    
    xRec = []; yRec = [];
    
else
    disp("Picking up from previous results");
    intRes = intmRslts.intRes;
    lbFval = intmRslts.lbFval;
    ubFval = intmRslts.ubFval;
    
    treeNode = intmRslts.treeNode;
    nodeNum = intmRslts.nodeNum;
    numNodes = intmRslts.numNodes;
    
    K = intmRslts.K;
    
    xRec = intmRslts.xRec; yRec = intmRslts.yRec;
    
    fvalCand = intmRslts.fvalCand;
    yValCand = intmRslts.yValCand;
    nodeNumCand = intmRslts.nodeNumCand;
    
end
ubLocal = inf;
ubFvalK = inf;
fprintf(header);
while( abs( lbFval - ubFvalK) > optimzOpts.tolGap && ...
        K < optimzOpts.iterLimt)
    % NodeNum - node number of smallest lower bound value -> parent node
    % for the new nodes
    
    % Level of current node (first level to backtrack for qualifying constraints)
    tLevel = treeNode(nodeNum).nodeLvl;
    % nLvl: Level for new node in the branch and bound tree
    nLvl = tLevel + 1;
    if(K == 1)
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
    end
    % Record the x and yK values during optimization process
    xRec = [xRec xPrim];
    yRec = [yRec [yK; lbFval]];
    % Do local search using updated x and yK values as initial values
    if(optimzOpts.lcaSearch)
        if(optimzOpts.xOption == 1)
            x0 = [xPrim(1:end-1); yK; xPrim(end)];
            lb = [x_lb(1:end-1); y_lb(1:end-1); x_lb(end)];
            ub = [x_ub(1:end-1); y_ub(1:end-1); x_ub(end)];
        else
            x0 = [yK; xPrim(1:end)];
            lb = [y_lb(1:end-1); x_lb(1:end)];
            ub = [y_ub(1:end-1); x_ub(1:end)];
        end
        [localTemp, fvalTemp, localFlag] = fmincon(localOptm.lcaObj,...
            x0, [], [], [], [], lb , ub, localOptm.lcaCons, nlinOptions);
        % fmincon finds better local optimal value
        if(localFlag > 0 && ubLocal > fvalTemp)
            ubLocal = fvalTemp;
            xLca = localTemp;
        end
        
        %Todo - how useful is following code
        % if better upper bound is found, solve the prime again to find
        % better Lagrange function
        %         if(abs(fval_ub - fval_ub_local) > 1e-8)
        %             yK = x_local(n_X + 1 : end-1);
        %             % Construct primal problem with new fixed yK value
        %             ineq = Inpt_Strct.constraints(1:end - (2 * n_X + 1));
        %             [A_prim, b_prim] = Vec_To_Mtx(yK,ineq,n_X + 1,n_Y);
        %             c = [zeros(n_X,1); 1];
        %             % Solve primal optimization problem
        %             if(lin_tolbx == 1)
        %                 [x, fval_prim, ~, ~, lambda] = cplexlp(c, A_prim, b_prim, [], [],...
        %                     x_lb, x_ub, [], options);
        %             else
        %                 [x, fval_prim, ~, ~, lambda] = linprog(c, A_prim, b_prim, [], [],...
        %                     x_lb, x_ub, [], options);
        %             end
        %             x_rec(:,K) = x;
        %             y_rec(:,K) = [yK; fval_lb];
        %         end
        
        % Gap between lcoal minimum and lower bound is smallet than
        % toletance value; exit the algtorithm with x, yK values
        if(abs(ubLocal - lbFval) < optimzOpts.tolGap)
            if(optimzOpts.xOption == 1)
                optmRslts.x = [xLca(1:numX); xLca(end)];
                optmRslts.yK = xLca(numX + 1 : end-1);
            else
                optmRslts.x = xLca(numX + 1 : end-1);
                optmRslts.yK = [xLca(1:numX); xLca(end)];
            end
            optmRslts.nodeNum = nodeNum;
            optmRslts.numNodes = numNodes;
            optmRslts.K =  K;
            optmRslts.xRec = xRec;
            optmRslts.yRec = yRec;
            optmRslts.treeNode = treeNode;
            optmRslts.lbFval = lbFval;
            optmRslts.ubFval = ubLocal;
            optmRslts.intRes = intRes;
            return
        end
        
    end
 
    % Update upper bound of the objective function
    ubFval = min([ubFval ubFvalK ubLocal]);
    
    
    % Lagrangian function substituted by xPrim
    % Example:  x2K is delta (not included in numX), and x2 is a
    % non-connected variable with Y.
    %      muK * (35.0 * x1K - 1.0 * x2K - 25.0 * y1 - 15.0 * x1K * y1 + 58.0)
    % Lag: 1 x (numY + 1). The entries are
    %      1 ~ numY: coefficients in front of Y_j
    %                Example: muK * (- 25.0 - 15.0 * x1K)
    %      (numY + 1): 
    %                Example: muK * (35.0 * x1K + 58.0)
    Lag = zeros(1,numY + 1);
    
    for i = 1 : numX
        Lag = Lag + xPrim(i) *  lambda.ineqlin' * coeffX(:,:,i);
    end
    Lag = Lag +  lambda.ineqlin' * constX;
    
    yLag = [Lag(1 : end - 1) -1];
    bLag = Lag(end);
    
    yGrad = zeros(numX, numY);
    bGrad = zeros(numX,1);
    
    for i = 1 : numX
        yGrad(i,:)  =  lambda.ineqlin' * coeffX(:,1 : end - 1, i);
        bGrad(i) = lambda.ineqlin' * coeffX(:, end, i) - lambda.lower(i) + lambda.upper(i);
    end
    
    delIdx = [];
    for i = 1 : numX
        if(max(abs(yGrad(i,:))) < 1e-5)
            delIdx = [delIdx i];
        end
    end
    
    conctIdx = setdiff( 1 : numX, delIdx);
    numConctX = length(conctIdx);
    yGrad = yGrad(conctIdx,:);
    bGrad = bGrad(conctIdx);
    
    %% Find qualifying constraints for relaxed subproblem
    AQual = []; bQual = [];
    % Current node is the first parent node of the new node
    pNode = nodeNum;
    while(tLevel > 1)
        AQual = [AQual; treeNode(pNode).ADualK];
        bQual = [bQual; treeNode(pNode).bDualK];
        tLevel = tLevel - 1;
        pNode = treeNode(pNode).pNode;
    end
    %% Binary tree for relaxed dual subproblem
    resIter = [];
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
            treePath = findpath(rexDulTree,1,botmNode(i));
            
            ADualK = zeros(1 + numConctX, numY + 1);
            bDualK = zeros(1 + numConctX, 1);
            
            xBl = zeros(numConctX, 1);
            % find value of Conncected X variables
            % find relaxed Lagrangian function and companying constraints
            % find value of Conncected X variables
            for j = 1 : numConctX
                if( rem(treePath(1 + j ),2) == 0 )
                    xBl(j) = x_ub(conctIdx(j));
                    ADualK(1 + j, :) = [yGrad(j,:) 0];
                    bDualK(1 + j) = -bGrad(j);
                else
                    % gradeint > 0 -> Ay + b > 0 -> -Ay < b
                    xBl(j) = x_lb(conctIdx(j));
                    ADualK(1 + j, :) = [-yGrad(j,:) 0];
                    bDualK(1 + j) = bGrad(j);
                end
            end
            % substitute x_B into Lagrange function to form relaxed dual problem
            ADualK(1,:) = yLag + [(xBl - xPrim(conctIdx))'  * yGrad 0];
            bDualK(1,:) = -(bLag + (xBl - xPrim(conctIdx))' * bGrad);
            
            if(nLvl == 2)
                ADual = ADualK;
                bDual = bDualK;
            else
                ADual = [AQual;ADualK];
                bDual = [bQual;bDualK];
            end
            % minimize mu_B
            c = [zeros(numY,1);
                1];
            
            if(optimzOpts.linTolbx == 1)
                [yTemp, ~, exitflag] = cplexlp(c, ADual, bDual, [], [],...
                    y_lb, y_ub, [], linOptions);
            else
                [yTemp, ~, exitflag] = linprog(c, ADual, bDual, [], [],...
                    y_lb, y_ub, [], linOptions);
            end
            
            if(exitflag > 0)
                yTemp(~any(ADual,1)) = yK(~any(ADual,1));
                % find zero column in A matrix -> correspondng yK value does
                % not update -> use value from last iteration
                yTemp(~any(ADual,1)) = yK(~any(ADual,1));
                % save yK and lower bound value
                dupSub = 1;
                %TODO: how to handle two subproblems with exactly same
                %answer so as to reduce computational time
                for j = 1 : size(resIter,2)
                    if(norm(yTemp - resIter(1 : end-1, j)) < 1e-7)
                        dupSub = 0;
                        CkNum = resIter(end,j);
                        
                        ref = [treeNode(CkNum).ADualK treeNode(CkNum).bDualK];
                        chk = [ADualK bDualK];
                        [~, dupIdx] = ismembertol(chk,ref,1e-7,'ByRows',true);
                        disIdx = setdiff(1 : size(ref,1),dupIdx);
                        del_idx = [];
                        for k = 1 : length(disIdx)
                            if(treeNode(CkNum).ADualK(disIdx(k),end) ~= -1)
                                del_idx = [del_idx disIdx(i)];
                            end
                        end
                        treeNode(CkNum).ADualK(del_idx,:) = [];
                        treeNode(CkNum).bDualK(del_idx) = [];
                        if(~ismembertol(chk(1,:),ref,1e-7,'ByRows',true))
                            treeNode(CkNum).ADualK = [treeNode(CkNum).ADualK; ADualK(1,:)];
                            treeNode(CkNum).bDualK = [treeNode(CkNum).bDualK; bDualK(1)];
                        end
                        break;
                    end
                end
                % make sure yK value is updated from the iteration
                if(norm(yTemp(1:numY) - yK) > 1e-7 && dupSub == 1)
                    
                    numNodes = numNodes + 1;
                    resIter = [resIter [yTemp; numNodes]];
                    treeNode(numNodes).nodeLvl = nLvl;
                    treeNode(numNodes).pNode = nodeNum;
                    treeNode(numNodes).yVal = yTemp(1:numY);
                    treeNode(numNodes).fval = yTemp(end);
                    % save ADualK and bDualK of new node
                    treeNode(numNodes).ADualK = ADualK;
                    treeNode(numNodes).bDualK = bDualK;
                end
            end
        end
    else
        % No conncected X variables
        numNodes = numNodes + 1;
        treeNode(numNodes).nodeLvl = nLvl;
        treeNode(numNodes).pNode = nodeNum;
        
        ADualK = yLag;
        bDualK = -bLag;
        
        treeNode(numNodes).ADualK = ADualK;
        treeNode(numNodes).bDualK = bDualK;
        
        if(nLvl == 2)
            ADual = ADualK;
            bDual = bDualK;
        else
            ADual = [AQual;ADualK];
            bDual = [bQual;bDualK];
        end
        
        c = [zeros(numY,1);
            1];
        
        if(optimzOpts.linTolbx == 1)
            [yTemp, ~, exitflag] = cplexlp(c, ADual, bDual, [], [],...
                y_lb, y_ub, [], linOptions);
        else
            [yTemp, ~, exitflag] = linprog(c, ADual, bDual, [], [],...
                y_lb, y_ub, [], linOptions);
        end
        
        if(exitflag > 0)
            yTemp(~any(ADual,1)) = yK(~any(ADual,1));
            treeNode(numNodes).yVal = yTemp(1:numY);
            treeNode(numNodes).fval = yTemp(end);
            % make sure yK value is updated from the iteration
            if(norm(yTemp(1:numY) - yK) > 1e-5)
                intRes = [intRes [yTemp; numNodes]];
            end
        end
    end
    
    intRes = [intRes resIter];
    % Fathom node with lb value larger than ub
    fathIdx = intRes(numY + 1,:) > (min(ubFval,ubLocal) + optimzOpts.tolGap);
    intRes(:,fathIdx) = [];
    
    % Find new lower bound value
    lbFval = min(intRes(numY + 1,:));
    
    % Find candidate with same lower bound value
    lbFvalIdx = find(abs(intRes(numY + 1,:) - lbFval) < 1e-7);
    
    nodeNum = intRes(end, lbFvalIdx(end));
    % Remove next yK value from all result group
    intRes(:,lbFvalIdx(end)) = [];
    
    fprintf(formatstr,K,min(ubFval,ubLocal),lbFval,ubFvalK,ubFvalK - lbFval,nodeNum,treeNode(nodeNum).pNode );
    
    
    if(K > 1)
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
    end
    K = K + 1;
    stem(1:numX,xPrim(1:numX),'o')
    ylim([min(x_lb(1:numX));max(x_ub(1:numX))]);
    xticks(1:1:numX);
    pause(0.5);
end
optmRslts.x = xPrim;
optmRslts.yK = yK;
optmRslts.nodeNum = nodeNum;
optmRslts.numNodes = numNodes;
optmRslts.K =  K;
optmRslts.xRec = xRec;
optmRslts.yRec = yRec;
optmRslts.treeNode = treeNode;
optmRslts.lbFval = lbFval;
optmRslts.ubFval = ubFval;
optmRslts.intRes = intRes;
end

