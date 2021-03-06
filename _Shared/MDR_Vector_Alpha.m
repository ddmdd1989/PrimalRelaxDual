function resPolyConstr = MDR_Vector_Alpha(K0, M0, K_j, lambdaExp, psiExp_m)

n_modes = length(lambdaExp);
n_alpha = size(K_j,3);
N = size(K0,1);
num_measDOFs = size(psiExp_m,1);
num_unmeasDOFs = N - num_measDOFs;

X = sym('x',[n_alpha + 1,1]);
Y = sym('y',[num_unmeasDOFs * n_modes,1]);

alpha = X(1 : n_alpha);
psiSim_u = reshape(Y, num_unmeasDOFs, n_modes);

psi = [psiExp_m; psiSim_u];
	
K = K0;

for i = 1:n_alpha
    K = K + alpha(i) * K_j(:,:,i);
end

for i = 1:n_modes
   resPoly((i-1) * N + 1 : i * N,1) = (K - lambdaExp(i) * M0) * psi(:,i);
end

n_resd = N * n_modes;

resPolyConstr = sym( zeros( 2 * n_resd, 1));

%% Modal dynamic residual constraints
for i = 1 : n_resd
    resPolyConstr(i,1) = resPoly(i) - X(end);
    resPolyConstr(n_resd + i,1) = -X(end) - resPoly(i);
end

