function [ineqPoly,bnd] = MDR_Vector_Psi(K0, M0, K_j, lambdaExp, psiExp_m, lb, ub)

n_modes = length(lambdaExp);
n_alpha = size(K_j,3);

N = size(K0, 1);
num_measDOFs = size(psiExp_m, 1);
num_unmeasDOFs = N - num_measDOFs;

%% Symbolic
X = sym('x',[num_unmeasDOFs * n_modes + 1,1]);
Y = sym('y',[n_alpha,1]);
numX = length(X);

alpha = Y(1:n_alpha);
psiSim_u = reshape(X(1:end-1),num_unmeasDOFs,n_modes);

psi = [psiExp_m;psiSim_u];

K = K0;

for i = 1:n_alpha
    
    K = K + alpha(i) * K_j(:,:,i);

end

for i = 1:n_modes

    polyResd((i-1) * N + 1 : i * N,1)= (K - lambdaExp(i) * M0) * psi(:,i);

end

n_eq = N * n_modes;
ineqPoly = sym(zeros(2 * n_eq,1));

for i = 1:n_eq
    
    ineqPoly(i,1) = polyResd(i) - X(end);
    ineqPoly(n_eq + i,1) = -X(end) - polyResd(i);
    
end

for i = 1:numX - 1
   
    bnd(i ,1) = sym(lb(i)) - X(i);
    bnd(numX +  i,1) = X(i) - sym(ub(i));

end
bnd(numX ,1) = -X(num_unmeasDOFs * n_modes + 1);




