
function generateplots(A, Ap, uf, u, ur, gtol, espai, alpha, beta, maxit)

n = length(A);
b = ones(n,1);

fprintf('Running GMRES-IR\n');
gmresir3(full(Ap),b,uf,u,ur,maxit,gtol);

fprintf('Running SGMRES-IR\n');
sgmresir3(full(Ap),b,uf,u,ur,maxit,gtol);

fprintf('Running GMRESIR_SPAI\n');
gmresir3_spai(full(A), b, espai, alpha, beta, uf, u, ur, maxit, gtol);

fprintf('Running SGMRESIR_SPAI\n');
sgmresir3_spai(full(A), b, espai, alpha, beta, uf, u, ur, maxit, gtol);

fprintf('Running GMRESIR_NP\n');
gmresir3_np(full(A), b, uf, u, ur, maxit, gtol);
