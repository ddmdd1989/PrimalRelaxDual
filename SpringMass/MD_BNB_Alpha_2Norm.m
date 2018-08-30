% clc;
clear;
% close all
warning('off');

%% Basis parameter

N = 4;                               % DOF of the whole structure
masses = 6 * ones(N, 1);             % kg
iniSpring = 35 * ones(N, 1);       % N/m
n_modes = 4;                        % number of measured mode shapes
dmgLoc = [1 3 ];
alpha_act = [0.2 -0.5];

tolGap = 1e-5;

actSpring = iniSpring;
for i = 1:length(dmgLoc)
    actSpring(dmgLoc(i)) = actSpring(dmgLoc(i)) * (1 + alpha_act(i));
end

%% Damaged strucuture
measDOFs = [1;2;3];
unmeasDOFs = 4;


M0 = makeM(masses, N);
K_act = makeK(actSpring, N);
[psiExp,lambdaExp] = eig(K_act, M0) ;
[lambdaExp,dummyInd] = sort((diag(lambdaExp)),'ascend') ;
lambdaExp = lambdaExp(1:n_modes);

psiExp = psiExp(:,dummyInd(1:n_modes));
psiExp_m = psiExp(measDOFs,:);
psiExp_u = psiExp(unmeasDOFs,:);

for i = 1:n_modes
    [~,I] = max(abs(psiExp_m(:,i)));
    psiExp_u(:,i) = psiExp_u(:,i) / psiExp_m(I,i);
    psiExp_m(:,i) = psiExp_m(:,i) /psiExp_m(I,i);
end

%% Influence Matrix
i = 1 ;
K_j(:,:,i) = zeros(N) ;
K_j(i,i,i) = iniSpring(i);

for i = 2 : N
    K_j(:,:,i) = zeros(N) ;
    K_j(i-1,i-1,i) = iniSpring(i);
    K_j(i-1,i,i) = -iniSpring(i);
    K_j(i,i-1,i) = -iniSpring(i);
    K_j(i,i,i) = iniSpring(i);
end

K_j = K_j(:,:,dmgLoc);


%% initial structure
M0 = makeM(masses, N);
K0 = makeK(iniSpring, N);
rordIdx = [measDOFs; unmeasDOFs];
K0 = K0(rordIdx, rordIdx);
M0 = M0(rordIdx, rordIdx);
K_j = K_j(rordIdx, rordIdx,:);


[psiSim,lambdaSim] = eig(K0,M0) ;

[lambdaSim,dummyInd] = sort((diag(lambdaSim)),'ascend') ;

psiSim_m = psiSim(1 : length(measDOFs),1:n_modes) ;
psiSim_u = psiSim(length(measDOFs) + 1 : end,1:n_modes) ;


for i = 1:n_modes
    [~,I] = max(abs(psiSim_m(:,i)));
    psiSim_u(:,i) = psiSim_u(:,i) / psiSim_m(I,i);
    psiSim_m(:,i) = psiSim_m(:,i) /psiSim_m(I,i);
end

n_alpha = size(K_j,3);
weight = ones(n_modes,1);
%% Bicvx variable
% alpha
X = sym('x',[n_alpha, 1]);
numX = length(X);  %% actual number of X variable

% Unmeasured entry of eigenvectors
Y = sym('y',[length(unmeasDOFs) * n_modes + 1,1]);
numY = length(Y) - 1;

x_lb = -ones(n_alpha,1);
x_ub =  ones(n_alpha,1);

y_lb =  [-3 * ones(length(unmeasDOFs) * n_modes,1); -inf];
y_ub =  [ 3 * ones(length(unmeasDOFs) * n_modes,1);  inf];


[ineq,bound] = MDR_Vector_2Norm(K0, M0, K_j, lambdaExp, psiExp_m, x_lb, x_ub);

optm = subs(ineq,[X;Y(1:end - 1)],[alpha_act';reshape(psiExp_u,length(unmeasDOFs) * n_modes,1)]);
y0 = zeros(length(unmeasDOFs) * n_modes,1);

prim = matlabFunction(subs(ineq,Y(1 : end - 1),y0),'Optimize',false,'Vars',{X});

nlinOptions = optimoptions('fmincon','Display','off',...
    'algorithm','interior-point','MaxFunEvals',50000,'MaxIter',10000);
fun_orig = @(x) Obj_dynamicresidual(x,K0,M0,K_j,lambdaExp, psiExp_m,weight);
polynomial = [ineq,;bound];
inptStrct = struct('y0', y0, 'polynomials',polynomial,'numX',numX);
iter_limit = 1000;
optStrct = struct('tolGap',tolGap,'iterLimt',iter_limit,'x_lb',x_lb,...
                    'x_ub',x_ub,'y_lb',y_lb,'y_ub',y_ub,'lcaSearch',1,...
                    'xOption',1);
loclStrct = struct('lcaObj',fun_orig);


optmRslts = Prim_Relx_Dual_2Norm(inptStrct,optStrct,loclStrct);

% [xIter, fval, localFlag,~, lambda] = fmincon(prim, zeros(n_alpha,1), [], [], [], [], x_lb , x_ub,[],nlinOptions);
% 
% lagFunc = ineq + [lambda.lower;lambda.upper;]' * bound;
% 
% conctX = []; xIdx = [];
% 
% for i = 1 : numX
%     [cx,tx] = coeffs(lagFunc,X(i));
%     if(~isempty(tx))
%        if(tx(end) == 1)
%             loopX = length(cx) - 1;
%        else
%            loopX = length(cx);
%        end
%         for j = 1:loopX
%             [~,tx] = coeffs(cx(j),Y);
%             if(~(length(tx) == 1 && tx(end) == 1) && ~isempty(tx))
%                 conctX = [conctX X(i)];
%                 xIdx = [xIdx i];    
%                 break;
%             end
%         end
%     end
% end
% 
% numConctX = length(conctX);
%    
% 
% rexDulTree = tree('root');
% for i = 1 : numConctX
%     for j = 1 : (2^i)
%         if(i == 1)
%             [rexDulTree, l1(j)] = rexDulTree.addnode(1, 0);
%         else
%             eval(['[rexDulTree, l' num2str(i) '(j)] = rexDulTree.addnode(l' num2str(i-1) '( (ceil(j/2)) ),0);']);
%         end
%     end
% end
% eval(['botmNode = l' num2str(numConctX) ';']);
% 
% gradLagCons = sym(zeros(numConctX,1));
% gradLagFunc = sym(zeros(1,numConctX));
% i = 1;
% treePath = findpath(rexDulTree,1,botmNode(i));
% 
% iterX = zeros(numConctX,1);
% 
% for j = 1 : numConctX
%     % differentiate Lagrangian function w.r.t connected variable
%     gradLagFunc(1,j) = subs(diff(lagFunc,conctX(j)),X,xIter);
%     if( rem(treePath(1 + j ),2) == 0 )
%         gradLagCons(j,1) = subs(diff(lagFunc,conctX(j)),X,xIter);
%         iterX(j) = x_ub(xIdx(j));
%     else 
%         gradLagCons(j,1) = -subs(diff(lagFunc,conctX(j)),X,xIter);
%         iterX(j) = x_lb(xIdx(j));
%     end
%     
% end
%      
%  ceq = 0;
%  linLag = (subs(lagFunc,X,xIter) + gradLagFunc * (iterX - xIter(xIdx))) - Y(end) ;
%  rexDual = [linLag;gradLagCons];
%  accpCons = matlabFunction(rexDual,ceq,'Outputs',{'cineq','ceq'},'Optimize',false,'Vars',{Y});
%  objDual = matlabFunction(Y(end),'Vars',{Y});
% 
%  [localTemp, fvalTemp, localFlag] = fmincon(objDual,...
%      [y0;0], [], [], [], [],y_lb , y_ub, accpCons, nlinOptions);
 

