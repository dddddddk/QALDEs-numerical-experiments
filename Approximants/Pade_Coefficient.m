function d = Pade_Coefficient(p)
% This function return the symmetric Pade coefficient to the approximants
% of Exponential function. The order is given by p.
d = zeros(p+1, 1);
d(1) = 1;
for i = 2:p+1
    d(i) = d(i-1) * (p - i + 2) / ((i-1) * (2*p - i + 2));
end
end