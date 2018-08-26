function [c,ceq] = Obj_dynamicresidual_constraints_infnorm(x,K0,M0,K_j,lambdExp,psiExp_m,weight)


n_alpha = size(K_j,3);
n_modes = size(psiExp_m,2) ;

N = size(M0,1) ;
psi = zeros(N,n_modes);
num_measDOFs = size(psiExp_m,1);
num_unmeasDOFs = N - num_measDOFs;

K = K0;
for i = 1 : n_alpha
    K = K + x(i) * K_j(:,:,i) ;
end

for i = 1 : n_modes
    psi(1 : num_measDOFs,i) = psiExp_m(:,i) ;
    psi(num_measDOFs + 1 : end,i) = x(n_alpha + (i-1) * num_unmeasDOFs + 1 : n_alpha + i * num_unmeasDOFs)';
end

for j = 1 : n_modes
    y_pos((j-1) * N + 1 : j * N, 1) = (K - lambdExp(j) * M0) * psi(:,j) * weight(j) - x(end);
    y_neg((j-1) * N + 1 : j * N, 1) = -x(end) - (K - lambdExp(j) * M0) * psi(:,j) * weight(j);
end

c = [y_pos;y_neg];
ceq = [];
end
