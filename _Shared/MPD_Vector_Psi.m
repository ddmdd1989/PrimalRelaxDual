function [ineqPoly,bnd] = MPD_Vector_Psi(K0, M0, K_j, lambdaExp, psiExp_m, q, eigWeight, lb, ub)


n_modes = length(lambdaExp);
n_alpha = size(K_j,3);
N = size(K0,1);
num_measDOFs = size(psiExp_m,1);
num_unmeasDOFs = N - num_measDOFs;
num_measDOFs_r = num_measDOFs - 1;

%% symbolic
X = sym('x',[( N - 1) * n_modes + 1,1]);
numX = length(X);

Y = sym('y',[n_alpha + n_modes,1]);

% Stiffness parameter
alpha = Y(1:n_alpha);

% Analytical eigvalue
lambdaSim = Y(n_alpha + 1 : end);

% Analytical eigenvector at measured DOF
psiSim_mR = reshape( X(1 : num_measDOFs_r * n_modes),...
                   num_measDOFs_r, n_modes);
% Analytical eigenvector at unmeasured DOF              
psiSim_u = reshape( X( num_measDOFs_r * n_modes + 1:end - 1),...
                   num_unmeasDOFs,n_modes);
 
% Organized analytical eigenvector at all DOF
psiSim = sym(zeros(N, n_modes));              

% Experimental eigenvector except for Qi entry
psiExp_mR = zeros(num_measDOFs - 1, n_modes);

for i = 1:n_modes
    % experimental eigenvector except for Qi entry
    measDOF_r = setdiff(1 : num_measDOFs, q(i));
    psiExp_mR(:,i) = psiExp_m(measDOF_r,i);
    
    % Assemble analytical eigenvector
    psiSim(measDOF_r, i) = psiSim_mR(:, i);
    psiSim(q(i), i) = 1;
    psiSim(num_measDOFs + 1 : end, i) = psiSim_u(:, i);
end

%% Assemble objective functiion
%  1. Eigenvalue difference
%  2. Eigenvector difference
%  3. weighed eigenvalue equation
K = K0;
for i = 1:n_alpha
    K = K + alpha(i) * K_j(:,:,i);
end

%% Assemble objective functiion
%  1. Eigenvalue difference
%  2. Eigenvector difference
%  3. weighed eigenvalue equation

lambdaDiff = sym(zeros(n_modes, 1));
psiDiff = sym(zeros((num_measDOFs - 1) * n_modes,1));
eigCons = sym(zeros(N * n_modes,1));

for i = 1:n_modes
    lambdaDiff(i,1) = (lambdaExp(i) - lambdaSim(i)) / lambdaExp(i);
    psiDiff((i-1) * num_measDOFs_r + 1 : i * num_measDOFs_r,1) = psiExp_mR(:,i) -  psiSim_mR(:,i);
    eigCons((i - 1) * N  + 1: i * N,1) = eigWeight * (K - lambdaSim(i) * M0) * psiSim(:,i);
end

diffPoly = [lambdaDiff; psiDiff; eigCons];

n_eq = length(diffPoly);
ineqPoly = sym(zeros(2 * n_eq,1));
for i = 1:n_eq
    ineqPoly(i,1) = diffPoly(i) - X(end);
    ineqPoly(n_eq + i,1) = -X(end) - diffPoly(i) ;
end

bnd = sym(zeros(2 * (numX) - 1,1));

for i = 1:(length(X) - 1)
    bnd(i ,1) = sym(lb(i)) - X(i);
    bnd(numX +  i,1) = X(i) - sym(ub(i));

end
bnd(numX ,1) = -X(end);












