clc;
clear all
data=load('fs_183_1.mtx');
A = spconvert(data);
n=length(A);
b=ones(n,1);
x = gmresir3(full(A), zeros(n,1), ones(n,1), .5, 20, 20, 1, 2, 4, 15, 1e-6,5,'fs_183_1e05');
drawnow
