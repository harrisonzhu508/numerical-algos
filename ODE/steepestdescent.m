function x_k = steepestdescent(A, b, x_0, c_1)
%%initialisation 
x_k = x_0;
complexity = 0;
p_k = -(A*x_k - b);
%%steepest descent algorithmc
while norm(p_k) > c_1
    %reinitialise a_k, the step, to 1 each iteration
    a_k = 1/2;
    p_k = -(A*x_k - b); %%define the steepest descent direction for each iteration
    %%choose a_k, the step, in this list such that the first Wolfe condition is
    %%satisfied
    while 0.5*(transpose(x_k + a_k*p_k)*A*(x_k + a_k*p_k) - transpose(x_k)*A*x_k) - transpose(b)*(x_k + a_k*p_k - x_k) > c_1*a_k*transpose(A*x_k - b)*p_k
        a_k = a_k/2;
    end
    x_k = x_k + a_k*p_k;
    complexity = complexity + 1;
end
complexity
end
