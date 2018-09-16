function [cineq, ceq] = MPDInfnormConstraints (x, K0, M0, K_j, lambdaExp, psiExp_m, q, eigWeight)

n_modes = length(lambdaExp);
n_alpha = size(K_j,3);
N = size(K0,1);
num_measDOFs = size(psiExp_m, 1);
num_measDOFs_r = num_measDOFs - 1;
num_unmeasDOFs = N - num_measDOFs;

% Stiffness parameter
alpha = x(1:n_alpha);
% Analytical eigvalue
lambdaSim = x(n_alpha + 1:n_alpha + n_modes);

numVar = n_alpha + n_modes;

% Analytical eigenvector at measured DOF
psiSim_mR = reshape(x(numVar + 1 : numVar + num_measDOFs_r * n_modes),...
                      num_measDOFs_r,n_modes);
                  
numVar = numVar + num_measDOFs_r * n_modes;

% Analytical eigenvector at unmeasured DOF
psiSim_u = reshape(x(numVar + 1: numVar + num_unmeasDOFs * n_modes),...
                     num_unmeasDOFs,n_modes);

% Organized analytical eigenvector at all DOF
psi_Sim = zeros(N,n_modes);

% Experimental eigenvector except for Qi entry
psi_mR = zeros(num_measDOFs - 1, n_modes);

for i = 1:n_modes
    % Assemble analytical eigenvector
    measDOF_r = setdiff(1 : num_measDOFs, q(i));
    psi_Sim(measDOF_r, i) = psiSim_mR(:,i);
    psi_Sim(q(i), i) = 1;
    psi_Sim(num_measDOFs + 1 : end, i) = psiSim_u(:, i);
    
    % experimental eigenvector except for Qi entry
    psi_mR(:,i) = psiExp_m(measDOF_r,i);
    
end
K = K0;
for i = 1 : n_alpha
    K = K + alpha(i) * K_j(:,:,i);
end

%% Assemble objective functiion
%  1. Eigenvalue difference
%  2. Eigenvector difference
%  3. weighed eigenvalue equation


lambdaDiff = zeros(n_modes,1);
psiDiff = zeros((num_measDOFs - 1) * n_modes,1);
eigCons = zeros(N * n_modes,1);

for i = 1:n_modes
    lambdaDiff(i,1) = (lambdaExp(i) - lambdaSim(i)) / lambdaExp(i);
    psiDiff((i-1) * num_measDOFs_r + 1 : i * num_measDOFs_r,1) = psi_mR(:,i) -  psiSim_mR(:,i);
    eigCons((i - 1) * N  + 1: i * N,1) = eigWeight * (K - lambdaSim(i) * M0) * psi_Sim(:,i);
end
polyDiff = [lambdaDiff;psiDiff;eigCons];
y_pos = zeros(length(polyDiff),1);
y_neg = zeros(length(polyDiff),1);
for j = 1 : length(polyDiff)
    y_pos(j,1) = polyDiff(j) - x(end);
    y_neg(j,1) = -x(end) - polyDiff(j);
end

cineq = [y_pos;y_neg];
ceq = [];
end
