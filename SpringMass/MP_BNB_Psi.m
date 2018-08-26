% clc;
clear;
close all

%% Basis parameter

N = 4;   % DOF of the whole structure
masses = 6 * ones(N, 1);                % kg
iniSpring = 35 * ones(N, 1);          % N/m

n_modes = 2;  % number of measured mode shapes
dmgLoc = [1 2 3 4];
alpha_act = [0.2; -0.5; 0.3; -0.3];

actSpring = iniSpring;

for i = 1:length(dmgLoc)
    actSpring(dmgLoc(i)) = actSpring(dmgLoc(i)).* (1 + alpha_act(i));
end

%% Damaged strucuture
measDOFs = [1;2;3];
unmeasDOFs = 4;

M0 = makeM(masses, N);
K_act = makeK(actSpring, N);

[psiExp, lambdaExp] = eig(K_act, M0) ;
[lambdaExp,dummyInd] = sort((diag(lambdaExp)),'ascend') ;
lambdaExp = lambdaExp(1 : n_modes);
modeIndex = 1 : n_modes;
psiExp = psiExp(:,dummyInd(modeIndex));

psiExp_m = psiExp(measDOFs, :);
psiExp_u = psiExp(unmeasDOFs, :);


psiExp_mR = zeros(length(measDOFs)-1, 1);
measDOFs_R = zeros(length(measDOFs)-1, 1);
q = zeros(n_modes,1);


for i = 1 : n_modes
    [~,q(i)] = max(abs(psiExp_m(:,i)));
    
    psiExp(:,i) = psiExp(:,i) / psiExp_m(q(i),i);
    psiExp_u(:,i) = psiExp_u(:,i) / psiExp_m(q(i),i);
    psiExp_m(:,i) = psiExp_m(:,i) / psiExp_m(q(i),i);  
    measDOFs_R(:,i) = setdiff(1 : length(measDOFs), q(i));
    psiExp_mR(:,i) = psiExp_m(measDOFs_R(:,i),i);
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
n_alpha = size(K_j,3);

alpha_ub = ones(n_alpha,1);
alpha_lb = -ones(n_alpha,1);


%% Initial Structure
M0 = makeM(masses, N);
K0 = makeK(iniSpring,N);

rordIdx = [measDOFs; unmeasDOFs;];
M0 = M0(rordIdx, rordIdx);
K0 = K0(rordIdx, rordIdx);
K_j = K_j(rordIdx, rordIdx,:);

[psiSim,lambdaSim] = eig(K0,M0);
[lambdaSim,dummyInd] = sort((diag(lambdaSim)),'ascend') ;
lambdaSim = lambdaSim(1:n_modes);
psiSim = psiSim(:,dummyInd(1:n_modes));
psiSim_m = psiSim(measDOFs,:) ;
psiSim_u = psiSim(unmeasDOFs,:);
psiSim_mR = zeros(length(measDOFs)-1,n_modes);
for i = 1:n_modes
    
    psiSim_u(:,i) = psiSim_u(:,i) / psiSim_m(q(i),i);
    psiSim_m(:,i) = psiSim_m(:,i) / psiSim_m(q(i),i);
   
    psiSim_mR(:,i) = psiSim_m(measDOFs_R(:,i),i);
    
end

%% Eigenvector except for qi-th entry + delta
X = sym('x',[(N-1) * n_modes + 1,1]);

x_lb = [-2 * ones((N - 1) * n_modes,1);   0];
x_ub = [ 2 * ones((N - 1) * n_modes,1); inf];
numX = (N - 1) * n_modes;

%% alpha + eigenvalue
Y = sym('y',[n_alpha + n_modes, 1]);
y_lb = [alpha_lb;lambdaExp * 0.5;-inf];
y_ub = [alpha_ub;lambdaExp * 1.5;inf];

eigWeight = 1 / norm(K0);

[ineqal,bound] = MPD_Vector_Psi(K0, M0, K_j, lambdaExp, psiExp_m, q, eigWeight, x_lb, x_ub);

Opt = [reshape(psiExp_mR,(length(measDOFs)- 1) * n_modes,1);reshape(psiExp_u,length(unmeasDOFs) * n_modes,1);0;alpha_act;lambdaExp;];

opt_ineqal = subs(ineqal,[X;Y],Opt);


polynomial = [ineqal;bound];

y0 = [zeros(n_alpha,1); lambdaSim;];

weight = ones(n_modes,1);

fun_orig = @(x)objfun_infnorm(x);
nonlcon = @(x)Obj_propertydiff_constraints_infnorm(x,K0,M0,K_j,lambdaExp,psiExp_m,q,eigWeight);



iter_limit = 1e5;
tolGap = 1e-5;
Inpt_Strct = struct('y0', y0, 'polynomials',polynomial,'numX',numX);
Opt_Strct = struct('tolGap',tolGap,'iterLimt',iter_limit,'x_lb',x_lb,...
                    'x_ub',x_ub,'y_lb',y_lb,'y_ub',y_ub,'lcaSearch',0,...
                    'linTolbx',2,'xOption',2);
Loc_Strct = struct('lcaObj',fun_orig,'lcaCons',nonlcon);



tic

final_result = Prim_Relx_Dual(Inpt_Strct,Opt_Strct,Loc_Strct);

t_plex = toc;


