function [besto, Q] = findbestreordering(A)

n = size(A,1);

% Find best reordering
[L,U] = lu(A);
nnzLU = nnz(L+U);

q = colperm(A);
[L,U] = lu(A(q,q));
nnzCP = nnz(L+U);

d = symrcm(A);
[L,U] = lu(A(d,d));
nnzRCM = nnz(L+U);

r = amd(A);
[L,U] = lu(A(r,r));
nnzMD = nnz(L+U);

p = dissect(A);
[L,U] = lu(A(p,p));
nnzND = nnz(L+U);

[~,ind] = min([nnzLU, nnzCP, nnzRCM, nnzMD, nnzND]);

switch ind
    case 1
        besto = 'nat';
        prm = [1:n]';
    case 2
        besto = 'clp';
        prm = q;
    case 3
        besto = 'rcm';
        prm = d;
    case 4
        besto = 'amd';
        prm = r;
    case 5
        besto = 'nds';
        prm = p;
end

I = eye(n);
Q = I(prm,:);