function ineqPoly = MPD_Vector_Alpha(K0, M0, K_j, lambdaExp, psi_m, q, eig_weight)

%x(1:n_alpha): stiffeness parameter (alpha);
%x(n_alpha + 1:n_alpha + numModes): Lambda_a;
%x(end): upper limit of one norm of modal property difference vector
%  
% Output:
%   ineqPoly
%   bndPoly

n_modes = length(lambdaExp);
n_alpha = size(K_j,3);
N = size(K0,1);
num_measDOF = size(psi_m,1);
num_unmeasDOF = N - num_measDOF;
num_measDOF_r = num_measDOF - 1;


X = sym('x',[n_alpha + n_modes + 1,1]);
n_X = length(X);
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
    eigCons((i - 1) * N  + 1: i * N,1) = eig_weight *  (K0 - lambdaSim(i) * M0) * psiSim(:,i);
end

diffPoly = [lambdaDiff; psiDiff; eigCons];

n_eq = length(diffPoly);
ineqPoly = sym(zeros(2 * n_eq,1));

for i = 1:n_eq
    ineqPoly(i,1) = diffPoly(i) - X(end);
    ineqPoly(n_eq + i,1) = -X(end) - diffPoly(i) ;
end












