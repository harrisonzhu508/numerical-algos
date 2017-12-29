function ps7_q3(N,tol)
format long
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
counter=0;
r=-b+A*x;
p=-r;
r2=1;

while (counter<N*N*2 && r2 >tol)%counter<=N*N
    
    r2=r'*r;
    alpha= r2/(p'*A*p);
    x= x + alpha* p;        %WATCH OUT!!!!
    r= r + alpha* A*p;
    beta= r'*r/r2;
    p = -r +beta*p;
    
    counter = counter +1;
end

q=1/(n+1):1/(n+1):1-1/(n+1);
true_ans=A\b;
x_final=x;
counter
r2
x(n/2)
plot(q,x_final,'DisplayName','x_{final}');hold on;plot(q,true_ans,'DisplayName','x_{true}');legend show;grid on
end
