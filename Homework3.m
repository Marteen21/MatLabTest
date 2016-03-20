%Ohne dich....
X0=93122186;
Y0=12345678;
X1=90091456;
Y1=51636639;
X2=91101022;
Y2=34626431;
P=99999989[

V0=[1,X0,mod(X0*X0,P)];
V1=[1,X1,mod(X1*X1,P)];
V2=[1,X2,mod(X2*X2,P)];

A=[V0;V1;V2]
B=[Y0;Y1;Y2]


C=A*transpose(A)
D=mod(inv(C)*A*B,P)