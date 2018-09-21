function ineqPoly = MPD_Vector_Alpha(K0, M0, K_j, lambdaExp, psi_m, q, eigWeight)

% function ineqPoly = MPD_Vector_Alpha(K0, M0, K_j, lambdaExp, psi_m, q, eigWeight)
%   (c) Yang Wang, Xinjun Dong (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2018
%
% This function using MATLAB symbolic toolbox to formualate the constraints
% of the eigenvector difference formulation for solving the problem with
% either Baron or P-RD optimization algorithm
% Input: 
%    K0 - nominal stiffness matrix (N x N), where N represents the number
%      of DOFs of the structure;
%    M0 - nominal mass matrix (N x N)
%    K_j - influence matrices for each updating variable (N x N x
%      n_alpha), where n_alpha represents the number of updating variables
%    lambdaExp - experimental eigenvalues (n_modes x 1), where n_modes
%      represents the number of experimental modes.
%    psi_m - experimental eigenvectors at measured DOFs (n_m x n_modes),
%      where n_m represents the number of measured DOFs
%    q - the entry with the maximum magnitude in psi_m (n_modes x 1)
%      In this formualtion, it is assumed that psi_m is normalized that the
%      maximum entry equa to 1, i.e. psi_m(q) =  1;
%    eigWeight - weighing factor for the eigenvalue equation (N * n_modes)
% Output:
%   ineqPoly - the constriants of the eigenvector difference formulation

n_modes = length(lambdaExp);
n_alpha = size(K_j,3);
N = size(K0,1);
num_measDOF = size(psi_m,1);
num_unmeasDOF = N - num_measDOF;
num_measDOF_r = num_measDOF - 1;


X = sym('x',[n_alpha + n_modes + 1,1]);
Y = sym('y',[(N-1) * n_modes,1]);

% Stiffness parameter
alpha = X(1:n_alpha);
% Analytical eigvalue
lambdaSim = X(n_alpha + 1 : n_alpha + n_modes);

% Analytical eigenvector at measured DOF
psiSim_mR = reshape(Y(1 : num_measDOF_r * n_modes),...
    num_measDOF_r,n_modes);

% Analytical eigenvector at unmeasured DOF
psiSim_u = reshape(Y(num_measDOF_r * n_modes + 1:end),...
    num_unmeasDOF,n_modes);

% Organized analytical eigenvector at all DOF
psiSim = sym(zeros(N,n_modes));

% Experimental eigenvector except for Qi entry
psi_mR = sym(zeros(num_measDOF - 1, n_modes));

for i = 1:n_modes
    measDOFs_r = setdiff(1 : num_measDOF, q(i));
    psi_mR(:,i) = psi_m(measDOFs_r,i);
    psiSim(measDOFs_r,i) = psiSim_mR(:,i);
    psiSim(q(i),i) = 1;
    psiSim(num_measDOF + 1 : end, i) = psiSim_u(:,i);
end

%% Assemble objective functiion
%  1. Eigenvalue difference
%  2. Eigenvector difference
%  3. weighed eigenvalue equation
for i = 1:n_alpha
    K0 = K0 + alpha(i) * K_j(:,:,i);
end

%% Assemble objective functiion
%  1. Eigenvalue difference
%  2. Eigenvector difference
%  3. weighed eigenvalue equation

lambdaDiff = sym(zeros(n_modes,1));
psiDiff = sym(zeros((num_measDOF - 1) * n_modes,1));
eigCons = sym(zeros(N * n_modes,1));

for i = 1:n_modes
    lambdaDiff(i,1) = (lambdaExp(i) - lambdaSim(i)) / lambdaExp(i);
    psiDiff((i-1) * num_measDOF_r + 1 : i * num_measDOF_r,1) = psi_mR(:,i) -  psiSim_mR(:,i);
    eigCons((i - 1) * N  + 1: i * N,1) = eigWeight(:,i) .*  (K0 - lambdaSim(i) * M0) * psiSim(:,i);
end

diffPoly = [lambdaDiff; psiDiff; eigCons];

n_eq = length(diffPoly);
ineqPoly = sym(zeros(2 * n_eq,1));

for i = 1:n_eq
    ineqPoly(i,1) = diffPoly(i) - X(end);
    ineqPoly(n_eq + i,1) = -X(end) - diffPoly(i) ;
end












