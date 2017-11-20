%script to create T_N matrix
function T_N = matrixTN(N)
Y = [2 -1; -1 2];
T_N = kron(eye(N/2,N/2),Y);
end
