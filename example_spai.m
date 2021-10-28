A=load('fs_183_1.m');
[m,n]=size(A);
for i=1:m
    j=A(i,1:2);
    we(j(1),j(2))=A(i,3);
end
A=we;
n=length(A);
b=ones(n,1);
%sir3(A,b,0.3,30,30,1,1,2,10);
gmresir3(A,b,0.1,17,17,1,1,2,10,1e-04);
