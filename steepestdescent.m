function x_k = steepestdescent(A, b, x_0, c_1, n, beta)
%%initialisation 
x_k = x_0;
%%define beta set and an empty set A
complexity = 0;
B = zeros(1,n);
for i = 1:n
    B(i) = beta^(i-1);
end
%%algorithm
for k = 1:n
    p_k = -(A*x_k - b);
    norm(p_k);
    if norm(p_k) <= c_1
        complexity
        return
    else
        %%choose a_k in this list such that the first Wolfe condition is
        %%satisfied
        Alpha = zeros(1,n);
        for i = 1:n
            Alpha(i)  =  ((0.5*transpose(x_k + B(i)*p_k)*A*(x_k + B(i)*p_k) - transpose(b)*(x_k + B(i)*p_k)) - (0.5*transpose(x_k)*A*x_k - transpose(b)*x_k) <= (c_1*B(i)*transpose(A*x_k - b)*p_k)) * B(i);
        end
        a_k = max(Alpha);
        x_k = x_k + a_k*p_k;
        complexity = complexity + 1;
    end
end
end
