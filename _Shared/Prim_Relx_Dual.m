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
%         Dimension 1: corresponding to each polynomial
%         Dimension 2: 1 -> numY - coeffecient of biliear terms
%                      numY + 1 - coeffecient of X linear term
%         Dimension 3: corresponding to each X variable
%       constX - terms not containing X variables (numPoly x (numY + 1))
%         Dimension 1: corresponding to each polynomial
%         Dimension 2: 1 -> numY - coeffecient of Y linear terms
%                      numY + 1 - constant term
%       coeffY - coeffecient of Y variables (numPoly x numX x numY)
%         Dimension 1: corresponding to each polynomial
%         Dimension 2: 1 -> numX - coeffecient of biliear terms excluding
%                      numX + 1 - coeffecient of Y linear term
%         Dimension 3: corresponding to each Y variable
%     constY - terms not containing Y variables (numPoly x (numX + 1))
%         Dimension 1: corresponding to each polynomial
%         Dimension 2: 1 -> numX - coeffecient of X linear terms
%                      numX + 1 - constant term
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
%   optmRslts - structure contains the final optimal result
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

x_lb = optimzOpts.x_lb; x_ub = optimzOpts.x_ub;
y_lb = optimzOpts.y_lb; y_ub = optimzOpts.y_ub;

numX = length(x_lb) - 1; % excluding the delta term
numY = length(optimzProb.y0);

coeffX = optimzProb.coeffX;
constX = optimzProb.constX;
coeffY = optimzProb.coeffY;
constY = optimzProb.constY;


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
        'MaxIterations', 1e8, 'ConstraintTolerance', 1e-5);
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
    %  -APoy : left side of qualifying Lagrangian and constraints
    %  -bPoy : right side of qualifying Lagrangian and constraints
    %    APoy * y <= bPoy;
    
    % root node in the branch and bound tree
    treeNode(1) = struct('nodeLvl', 1, 'pNode', 0, 'yVal', optimzProb.y0,...
        'APoy', [], 'bPoy', []);
    nodeNum = 1;
    numNode = 1;
    iter = 1;
    
    xRec = []; yRec = [];
    
else
    disp("Picking up from previous results");
    intRes = intmRslts.intRes;
    lbFval = intmRslts.lbFval;
    ubFval = intmRslts.ubFval;
    
    treeNode = intmRslts.treeNode;
    nodeNum = intmRslts.nodeNum;
    numNode = intmRslts.numNode;
    
    iter = intmRslts.iter;
    
    xRec = intmRslts.xRec; yRec = intmRslts.yRec;
    
    fvalCand = intmRslts.fvalCand;
    yValCand = intmRslts.yValCand;
    nodeNumCand = intmRslts.nodeNumCand;
    
end
ubLocal = inf;
fprintf(header);
ubIter = inf;
while( abs( lbFval - ubIter) > optimzOpts.tolGap && ...
        iter < optimzOpts.iterLimt)
    % NodeNum - node number of smallest lower bound value -> parent node
    % for the new nodes
    % Level for new node in the branch and bound tree
    % -> level of parent node + 1
    nLvl = treeNode(nodeNum).nodeLvl + 1;
    % Level of current node (first level to backtrack for qualifying constraints)
    tLevel = treeNode(nodeNum).nodeLvl;
    % Fixed y value for solving primal problem
    y = treeNode(nodeNum).yVal;
    % Construct primal problem with new fixed y value
    A_prim = zeros(size(coeffX,1),numX + 1);
    for i = 1: numX + 1
        A_prim(:,i) = coeffX(:,:,i) * [y; 1];
    end
    b_prim = -constX * [y; 1];
    c = [zeros(numX,1); 1];
    % Solve primal optimization problem
    if(optimzOpts.linTolbx == 1)
        [x, primFval, ~, ~, lambda] = cplexlp(c, A_prim, b_prim, [], [],...
            x_lb, x_ub, [], linOptions);
    else
        [x, primFval, ~, ~, lambda] = linprog(c, A_prim, b_prim, [], [],...
            x_lb, x_ub, [], linOptions);
    end
    % Record the x and y values during optimization process
    xRec = [xRec x];
    yRec = [yRec [y; lbFval]];
    % Update upper bound of the objective function
    ubFval = min([ubFval, primFval]);
    % Do local search using updated x and y values as initial values
    if(optimzOpts.lcaSearch)
        if(optimzOpts.xOption == 1)
            x0 = [x(1:end-1); y; x(end)];
            lb = [x_lb(1:end-1); y_lb(1:end-1); x_lb(end)];
            ub = [x_ub(1:end-1); y_ub(1:end-1); x_ub(end)];
        else
            x0 = [y; x(1:end)];
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
        %             y = x_local(n_X + 1 : end-1);
        %             % Construct primal problem with new fixed y value
        %             ineq = Inpt_Strct.constraints(1:end - (2 * n_X + 1));
        %             [A_prim, b_prim] = Vec_To_Mtx(y,ineq,n_X + 1,n_Y);
        %             c = [zeros(n_X,1); 1];
        %             % Solve primal optimization problem
        %             if(lin_tolbx == 1)
        %                 [x, fval_prim, ~, ~, lambda] = cplexlp(c, A_prim, b_prim, [], [],...
        %                     x_lb, x_ub, [], options);
        %             else
        %                 [x, fval_prim, ~, ~, lambda] = linprog(c, A_prim, b_prim, [], [],...
        %                     x_lb, x_ub, [], options);
        %             end
        %             x_rec(:,iter) = x;
        %             y_rec(:,iter) = [y; fval_lb];
        %         end
        
        % Gap between lcoal minimum and lower bound is smallet than
        % toletance value; exit the algtorithm with x, y values
        if(abs(ubLocal - lbFval) < optimzOpts.tolGap)
            if(optimzOpts.xOption == 1)
                optmRslts.x = [xLca(1:numX); xLca(end)];
                optmRslts.y = xLca(numX + 1 : end-1);
            else
                optmRslts.x = xLca(numX + 1 : end-1);
                optmRslts.y = [xLca(1:numX); xLca(end)];
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
    % Update upper bound of objective function value
    % ubFval = min(ubFval, ubLocal);
    % Upper bound of objective function value for current iteration
    ubIter = primFval;
    
     yLag = zeros(1,numY + 1);
    yLag(end) = -1;
    for i = 1 : numY
        yLag(1, i) = lambda.ineqlin' * coeffY(:,:,i) * [x; 1];
    end
    bLag = lambda.ineqlin' * constY(:,[1:end-2 end]) * [x(1:end - 1); 1]...
           + lambda.lower' * (x_lb - x) - lambda.upper(1 : numX)' * (x(1 : numX) - x_ub(1 : numX));

    yGrad = zeros(numX, numY);
    bGrad = zeros(numX,1);
    for i = 1 : numX 
        yGrad(i,:)  =  lambda.ineqlin' * coeffX(:,1 : end - 1, i);
        bGrad(i) = lambda.ineqlin' * coeffX(:, end, i) - lambda.lower(i) + lambda.upper(i);
    end
        
    delIdx = [];
    for i = 1 : numX
        if(sum(yGrad(i,:)) == 0)
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
        AQual = [AQual; treeNode(pNode).APoy];
        bQual = [bQual; treeNode(pNode).bPoy];
        tLevel = tLevel - 1;
        pNode = treeNode(pNode).pNode;
    end
    %% Binary tree for relaxed dual subproblem
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
            
            APoy = zeros(1 + numConctX, numY + 1);
            
            bPoy = zeros(1 + numConctX, 1);
            
            treePath = findpath(rexDulTree,1,botmNode(i));
            
            iterX = zeros(numConctX, 1);
            % find value of Conncected X variables
            % find relaxed Lagrangian function and companying constraints
           % find value of Conncected X variables
            for j = 1 : numConctX
                if( rem(treePath(1 + j ),2) == 0 )
                    iterX(j) = x_ub(conctIdx(j));
 
                    APoy(1 + j, :) = [yGrad(j,:) 0];
                    bPoy(1 + j) = -bGrad(j);
                else
                    % gradeint > 0 -> Ay + b > 0 -> -Ay < b
                    iterX(j) = x_lb(conctIdx(j));
                    APoy(1 + j, :) = [-yGrad(j,:) 0];
                    bPoy(1 + j) = bGrad(j);
                end
            end
            % substitute x_B into Lagrange function to form relaxed dual problem
                       
            APoy(1,:) = yLag + [(iterX - x(conctIdx))'  * yGrad 0];
            bPoy(1,:) = -(bLag + (iterX - x(conctIdx))' * bGrad);
                       
            % save APoy and bPoy of new node
            treeNode(numNode).APoy = APoy;
            treeNode(numNode).bPoy = bPoy;
            
            if(nLvl == 2)
                ADual = APoy;
                bDual = bPoy;
            else
                ADual = [AQual;APoy];
                bDual = [bQual;bPoy];
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
                % find zero column in A matrix -> correspondng y value does
                % not update -> use value from last iteration
                yTemp(~any(ADual,1)) = y(~any(ADual,1));
                % save y and lower bound value
                treeNode(numNode).yVal = yTemp(1:numY);
                treeNode(numNode).fval = yTemp(end);
                % make sure y value is updated from the iteration
                if(norm(yTemp(1:numY) - y) > 1e-7)
                    intRes = [intRes [yTemp; numNode]];
                end
            end
        end
    else
        % No conncected X variables
        numNode = numNode + 1;
        treeNode(numNode).nodeLvl = nLvl;
        treeNode(numNode).pNode = nodeNum;
       
        APoy = yLag;
        bPoy = -bLag;
       
        treeNode(numNode).APoy = APoy;
        treeNode(numNode).bPoy = bPoy;
        
        if(nLvl == 2)
            ADual = APoy;
            bDual = bPoy;
        else
            ADual = [AQual;APoy];
            bDual = [bQual;bPoy];
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
            yTemp(~any(ADual,1)) = y(~any(ADual,1));
            treeNode(numNode).yVal = yTemp(1:numY);
            treeNode(numNode).fval = yTemp(end);
            % make sure y value is updated from the iteration
            if(norm(yTemp(1:numY) - y) > 1e-5)
                intRes = [intRes [yTemp; numNode]];
            end
        end
    end
    
    % Fathom node with lb value larger than ub
    fathIdx = intRes(numY + 1,:) > (min(ubFval,ubLocal) + optimzOpts.tolGap);
    intRes(:,fathIdx) = [];
    
    % Find new lower bound value
    lbFvalIter = min(intRes(numY + 1,:));
    
    % Find candidate with same lower bound value
    lbFvalIdx = find(abs(intRes(numY + 1,:) - lbFvalIter) < 1e-7);
    
    % Pick y candidates from current iteration achieve the same lower
    % bound value
    yValIter = intRes(1 : numY, lbFvalIdx);
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
    c = [zeros(numX,1); 1];
    for i = 1:size(yValIter,2)
        
        ACand = zeros(size(coeffX,1),numX + 1);
        for j = 1: numX + 1
            ACand(:,j) = coeffX(:,:,j) * [yValIter(:,i); 1];
        end
        bCand = -constX * [yValIter(:,i); 1];

        if(optimzOpts.linTolbx == 1)
            [~, fvalCanIter(i)] = cplexlp(c, ACand, bCand, [], [],...
                x_lb, x_ub, [], linOptions);
        else
            [~, fvalCanIter(i)] = linprog(c, ACand, bCand, [], [],...
                x_lb, x_ub, [], linOptions);
        end
    end
    % Add upper bound value from new y candidate to candidate group
    fvalCand = [fvalCand fvalCanIter];
    yValCand = [yValCand yValIter];
    nodeNumCand = [nodeNumCand nodeNumIter];
    
    % Find smallest upper bound from all y (including current and previous
    % iterations) to be used for next iteration
    [~,min_idx] = min(fvalCand);
    %     [~,min_idx] = max(fvalCand);
    nodeNum = nodeNumCand(min_idx);
    
    % Remove next y value from candidate group
    fvalCand(min_idx) = [];
    yValCand(:,min_idx) = [];
    nodeNumCand(min_idx) = [];
    
    % Remove next y value from all result group
    intRes(:,intRes(end,:) == nodeNum) = [];
    
    fprintf(formatstr,iter,min(ubFval,ubLocal),lbFval,ubIter,ubIter - lbFval,nodeNum,treeNode(nodeNum).pNode );
    iter = iter + 1;
    stem(1:numX,x(1:numX),'o')
    ylim([min(x_lb(1:numX));max(x_ub(1:numX))]);
    xticks(1:1:numX);
    pause(0.5);
end
optmRslts.x = x;
optmRslts.y = y;
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

