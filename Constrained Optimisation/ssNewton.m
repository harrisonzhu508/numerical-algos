function [x,lambda]=ssNewton(N,iter,gamma)
A = (N+1)*(N+1)*(sparse(2:N,1:N-1,-1,N,N) + sparse(1:N,1:N,2,N,N) + sparse(1:N-1,2:N,-1,N,N));
mat1 = horzcat(A,-speye(N));
mat2 = horzcat(speye(N),sparse(N,N));
mat  = vertcat(mat1,mat2);
b = -ones([N,1]);
c = -0.05*ones([N,1]);
rhs = vertcat(b,c);
full(mat);
rhs;
lambda1=zeros([N,1]);
x1=zeros([N,1]);
xlambda1=vertcat(x1,lambda1);
for n=1:iter
  matmod=mat;
  rhsmod=rhs;
  for i=1:N
    if (xlambda1(N+i)<gamma*(xlambda1(i)-c(i)))
      matmod(i,N+i)=0;
      matmod(N+i,i)=0;
      matmod(N+i,N+i)=1;
      rhsmod(N+i)=0;
    end
    end
  xlambda2=matmod\rhsmod;
discrep=norm(xlambda2-xlambda1)/norm(xlambda2);
  xlambda1=xlambda2;
  fprintf ('Iter: %d Discrepancy: %f \n',n,discrep);
  if (discrep<0.0001) 
    break
  end
end
x=xlambda2(1:N);
lambda=xlambda2(N+1:2*N);


end