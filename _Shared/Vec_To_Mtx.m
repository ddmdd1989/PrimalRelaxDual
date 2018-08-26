function [A,b] = Vec_To_Mtx(y,polynomial_ineqal,n_X,n_Y)


X = sym('x',[n_X,1]);
Y = sym('y',[n_Y,1]);

polynomial_ineqal = subs(polynomial_ineqal,Y,y);

A = zeros(length(polynomial_ineqal),length(X));
b = zeros(length(polynomial_ineqal),1);


for i = 1:length(polynomial_ineqal)
    [cx,tx] = coeffs(polynomial_ineqal(i));
    cx = double(cx);
    [~,ia,ib] = intersect(tx,X);
    A(i,ib) = cx(ia);
    if(~isempty(find(tx == 1, 1)))
        b(i,1) = -cx(end);            
    end
end
end

