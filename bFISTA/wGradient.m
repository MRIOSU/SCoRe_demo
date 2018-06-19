function grad = wGradient(y,x,A,At)
% computes the gradient of the data fidelity
grad = 2*At(A(x) - y);