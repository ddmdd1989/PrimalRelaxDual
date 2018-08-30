function [polynomial,bnd] = MDR_Vector_2Norm(K0, M0, K_j, lambdaExp, psiExp_m, lb, ub)

n_modes = length(lambdaExp);
n_alpha = size(K_j,3);
N = size(K0,1);
num_measDOFs = size(psiExp_m,1);
num_unmeasDOFs = N - num_measDOFs;

X = sym('x',[n_alpha,1]);
Y = sym('y',[num_unmeasDOFs * n_modes,1]);

alpha = X(1 : n_alpha);
psiSim_u = reshape(Y, num_unmeasDOFs, n_modes);

psiMix = [psiExp_m; psiSim_u];
	
K = K0;

for i = 1:n_alpha
    K = K + alpha(i) * K_j(:,:,i);
end

for i = 1:n_modes
   resPoly((i-1) * N + 1 : i * N,1) = (K - lambdaExp(i) * M0) * psiMix(:,i);
end

polynomial = resPoly.' * resPoly; 

%% Boundary constraints
bnd = sym(zeros(n_alpha * 2, 1));
for i = 1 : n_alpha
    bnd(i, 1) = sym(lb(i)) - X(i);
    bnd(n_alpha + i, 1) = X(i) - sym(ub(i));
end








