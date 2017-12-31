function ps6_q3(N,iter)
A = (N+1)*(N+1)*(sparse(2:N,1:N-1,-1,N,N) + sparse(1:N,1:N,2,N,N) + sparse(1:N-1,2:N,-1,N,N));
b = 2*ones([N,1]);
x1=ones([N,1]);
matmod= A;
rhs=b;
for n=1:iter
  
  for i=1:N
    matmod(i,i) = A(i,i) +4*(x1(i))^3;
    rhs(i) = b (i) + 3*(x1(i))^4;
  end
  x2=matmod\rhs;
  discrep=norm(x2-x1)/norm(x2);
  x1=x2;
  fprintf ('Iter: %d Discrepancy: %f \n',n,discrep);
  if (discrep<0.0001) 
    break
  end
end

u=(0:1/(N-1):1);
x_newton=x2(1:N);


fun =@(z) 0.5*z'*A*z + 0.2*z.^5'*ones(N,1) -2*z'*ones(N,1);
[x_spgf] = fminunc(fun,ones([N,1]),opts);


plot(u,x_newton,'-b*',u,x_spgf,'-.r'),grid on


end
