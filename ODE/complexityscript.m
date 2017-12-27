function complexity = complexityscript(N)

A = sparse(((N+1)^2)*matrixTN(N));
b = ones(N,1);
x0 = zeros(N,1);
%%x1 = rand(N,1);
xstar = steepestdescent(A, b, x0, 10^(-4))
%%0.5*transpose(xstar)*A*xstar - transpose(b)*xstar
%%0.5*transpose(x1)*A*x1 - transpose(b)*x1

end