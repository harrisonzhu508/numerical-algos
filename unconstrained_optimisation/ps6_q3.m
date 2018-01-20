%we approximate -u''(x) + u(x)^4 = 2, 0 < x < 1 
%               u(0) = u(1) = 0
%via finite difference and then apply the Newton method to solve the linear
%system.
%for this problem, h = 1/(n+1) will be the time step and H will denote the
%Hessian of the function 
%----------------------------------------------------------------------------

function ps6_q3(N,iter)
A = (N+1)*(N+1)*(sparse(2:N,1:N-1,-1,N,N) + sparse(1:N,1:N,2,N,N) + sparse(1:N-1,2:N,-1,N,N));
b = 2*ones([N,1]);
x1=ones([N,1]);
hess=zeros(N,N);
gradf=zeros(N,1);
x1(1)=0;x1(N)=0;
for n=1:iter
  
    
  
  for i=1:N
      hess(i,i) =  4*(x1(i))^3;
      if i == 1
          gradf(i) = 2*x1(1) - x1(2) + (x1(i))^4 - b(i);
      elseif i == N
          gradf(i) = -x1(N-1) + 2*x1(N) + (x1(i))^4 - b(i);
      else   
          gradf(i) = -x1(i-1)-x1(i+1)+2*x1(i) + (x1(i))^4 - b(i);
      end
  end
  hess = A + hess;
  pk=-hess\gradf;
  x2= x1 + pk;
  discrep=norm(pk);
  x1=x2;
  fprintf ('Iter: %d Discrepancy: %f \n',n,discrep);
  if (discrep<0.0001) 
    break
  end
end
u=(0:1/(N-1):1);
x_newton=x2(1:N);

fun =@(z) 0.5*z'*A*z + 0.2*z.^5'*ones(N,1) -2*z'*ones(N,1);
[x_spgf] = fminunc(fun,ones([N,1]));


plot(u,x_newton,'-b*',u,x_spgf,'-.r'),grid on


end