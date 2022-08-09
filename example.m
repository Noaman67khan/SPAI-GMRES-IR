%Example script for comparing GMRES-IR, GMRESIR_SPAI and GMRESIR_NP (with 3 precisions)

clc;
clear all
load('pores_3');
A = Problem.A;

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
        Ap = A;
        besto = 'nat';
    case 2
        Ap = A(q,q);
        besto = 'clp';
    case 3
        Ap = A(d,d);
        besto = 'rcm';
    case 4
        Ap = A(r,r);
        besto = 'amd';
    case 5
        Ap = A(p,p);
        besto = 'nds';
end

besto
% Run experiments
generateplots(A, Ap, 1, 2, 4, 1e-8, .3, 5, 8, 15)

