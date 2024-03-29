
function generateplots(A, Q, uf, u, ur, gtol, espai, alpha, beta, maxit)

n = length(A);
b = ones(n,1);
b = b./norm(b);


% fprintf('Running GMRES-IR\n');
% gmresir3(full(Ap),b,uf,u,ur,maxit,gtol);

%('Running SGMRES-IR, (%d, %d, %d)\n', uf, u, ur);
%gmresir3(full(A),Q,b,uf,u,ur,maxit,gtol);

% fprintf('Running GMRESIR_SPAI\n');
% gmresir3_spai(full(A), b, espai, alpha, beta, uf, u, ur, maxit, gtol);
% 
fprintf('Running SGMRESIR_SPAI_scaling, (%d, %d, %d)\n', uf, u, ur);
sgmresir3_spai_scaling(full(A), b, espai, alpha, beta, uf, u, ur, maxit, gtol);

fprintf('Running SGMRESIR-SPAI, (%d, %d, %d)\n', uf, u, ur);
sgmresir3_spai(full(A), b, espai, alpha, beta, uf, u, ur, maxit, gtol);

fprintf('Running GMRESIR-NP, (%d, %d, %d)\n', uf, u, ur);
gmresir3_np(full(A), b, uf, u, ur, maxit, gtol);
