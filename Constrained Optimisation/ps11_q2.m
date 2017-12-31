function [u,x,lambda,s]=ps11_q2(N,iter,epsilon)

A = (N+1)*(N+1)*(sparse(2:N,1:N-1,-1,N,N) + sparse(1:N,1:N,2,N,N) + sparse(1:N-1,2:N,-1,N,N));
mat1 = horzcat(A,-speye(N),sparse(N,N));
mat2 = horzcat(speye(N),sparse(N,N),speye(N));
%mat3 = horzcat(sparse(N,N),sparse(diag(s1)),sparse(diag(lambda1)));
%mat  = vertcat(mat1,mat2,mat3);
b = -ones([N,1]);
c = -0.05*ones([N,1]);
%slambda = s1.*lambda1;
%rhs = vertcat(b,c,slambda);
%full(mat);
%rhs;
u=(0:1/(N-1):1);
x1=A\b;
lambda1=ones([N,1])*epsilon;
s1=ones([N,1])*epsilon;
%full(mat3)
%full(matmod)

xlambdas1=vertcat(x1,lambda1,s1);
%pause;
for n=1:iter
  
  for i=1:N
      if (xlambdas1(N+i)<epsilon)% && i>1)
          xlambdas1(N+i)=epsilon;
      end
      if (xlambdas1(2*N+i)<epsilon)% && i>1)
          xlambdas1(2*N+i)=epsilon;
      end
  end
  
  mat3 = horzcat(sparse(N,N),sparse(diag(xlambdas1(2*N+1:end))), ...
                          sparse(diag(xlambdas1(N+1  :2*N))) );
  slambda = xlambdas1(2*N+1:end).*xlambdas1(N+1  :2*N);
  
  
  matmod = vertcat(mat1,mat2,mat3);
  rhsmod = vertcat(b,c,slambda);
  full(matmod);
  full(rhsmod);
  %pause;
  
  xlambdas2=matmod\rhsmod;
  discrep=norm(xlambdas2(2*N+1:end)-xlambdas1(2*N+1:end))/norm(xlambdas2(2*N+1:end));
  xlambdas1=xlambdas2;
 % plot(u,xlambdas2(1:N)),grid on, hold on
  %pause;
  fprintf ('Iter: %d Discrepancy: %f \n',n,discrep);
  if (discrep<0.0001) 
    break
  end
end
%full(mat3)
%full(matmod)

x=xlambdas2(1:N);
lambda=xlambdas2(N+1:2*N);
s=xlambdas2(2*N+1:end);
plot(u,xlambdas2(1:N)),grid on


end
