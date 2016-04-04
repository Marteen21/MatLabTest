A = randi(7,6,5)
rref(A)
M=A(:,1:4)
b=A(:,5)
[Q,R] = qr(M,0)
x=inv(R)*Q'*b
x1=R\Q'*b
x2=pinv(M)*b