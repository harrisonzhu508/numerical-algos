function ps5_q3(N,tol)

T=2*eye(N);
for i=(2:N)
    T(i,i-1)=-1;
    T(i-1,i)=-1;
end


T=sparse(T);

%{
%a
n=N;    
A=(N+1)^2*T;
%}



%b
n=N*N;  
A=sparse(kron(eye(N),T))+sparse(kron(T,eye(N)));
A=(N+1)^2*sparse(A);


%{
%c
n=N*N*N;
A=kron(sparse(kron(eye(N),eye(N))) ,T)+ kron(sparse(kron(eye(N),T)),eye(N)) + kron(sparse( kron(T,eye(N))),eye(N));
A=(N+1)^2*sparse(A);
%}

b=ones(n,1);
x=ones(n,1);


%setting
c=0.0001;
beta=0.5;
counter=0;
p2=1;

while (counter<=N*N*2 && p2>tol)%NORM>tol)
    counter = counter +1;
    %f=ps5_q3_f(A,b,x);  %function values
    p=A*x-b;            %p = + grad(f)
    p2=p'*p;       %norm of p

    %{
    %Armijo rule
    alpha=1;        
    while ( ps5_q3_f(A,b,x-alpha*p) > f - c*alpha/NORM)
        alpha = alpha *beta;
    end
    %}
    
    alpha=p2/(p'*A*p);
    x = x - alpha*p;
end
q=1/(n+1):1/(n+1):1-1/(n+1);
true_ans=A\b;
x_final=x;
counter
p2
plot(q,x_final,'DisplayName','x_{final}');hold on;plot(q,true_ans,'DisplayName','x_{true}');legend show;grid on
end
