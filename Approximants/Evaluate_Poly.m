function val = Evaluate_Poly(x, d)
% d: the polynomial coefficients
% x: the input value, scalar or matrix
sz = size(x); p = length(d);
I = eye(sz);
val = d(p) * x + d(p-1) * I; 
for i = (p-2):-1:1
    val = x * val + d(i) * I; 
end
end