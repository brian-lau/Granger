l = 4;
m = 3;

A = randn(l,l);
B = randn(l,m);
C = randn(m,l);
D = randn(m,m);

A_hat = [A,B;C,D]

F2 = inv(D - C*inv(A)*B)
[inv(A)*(eye(l) + B*F2*C*inv(A)) , -inv(A)*B*F2 ; -F2*C*inv(A) , F2]

G2 = inv(A - B*inv(D)*C)
