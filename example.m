%Example script for comparing GMRES-IR, GMRESIR_SPAI and GMRESIR_NP (with 3 precisions)

% clc;
% clear all
load('steam3');
A = Problem.A;
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

besto
I = eye(n);
Q = I(prm,:);

% Run experiments
generateplots(A, Q, 1, 2, 4, 1e-8, .5, 10, 8, 15)

