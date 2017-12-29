function ps5_q3(N,tol)

T=2*eye(N);
for i=(2:N)
    T(i,i-1)=-1;
    T(i-1,i)=-1;
end

%{
%a
n=N;    
A=(N+1)^2*T;
%}

%{
%b
n=N*N;  
A=kron(eye(N),T)+kron(T,eye(N));
A=(N+1)^2*A;
%}


%c
n=N*N*N;
A=kron(kron(eye(N),eye(N)),T)+kron(kron(eye(N),T),eye(N))+kron(kron(T,eye(N)),eye(N));
A=(N+1)^2*A;


b=(ones(1,n))';
x=ones(n,1);

%setting
c=0.0001;
beta=0.5;
NORM=1;
counter=0;

while (counter<=N*N*2)%NORM>tol)
    counter = counter +1;
    f=ps5_q3_f(A,b,x);  %function values
    p=A*x-b;            %p = + grad(f)
    NORM=norm(p);       %norm of p

    %{
    %Armijo rule
    alpha=1;        
    while ( ps5_q3_f(A,b,x-alpha*p) > f - c*alpha/NORM)
        alpha = alpha *beta;
    end
    %}
    
    alpha=NORM*NORM/(p'*A*p);
    x = x - alpha*p;
end
true_ans=A\b;
x_final=x;
plot(x_final,'DisplayName','x_{final}');hold on;plot(true_ans,'DisplayName','x_{true}');legend show;grid on; hold off
end
