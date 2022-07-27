%Example script for comparing GMRES-IR, GMRESIR_SPAI and GMRESIR_NP (with 3 precisions)

clc;
clear all
data = load('orsreg_1');
A = spconvert(data.Problem.A);

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
    case 2
        Ap = A(q,q);
    case 3
        Ap = A(d,d);
    case 4
        Ap = A(r,r);
    case 5
        Ap = A(p,p);
end

% Run experiments
generateplots(A, Ap, 0, 1, 2, 1e-4, .5, 50, 50, 15)

