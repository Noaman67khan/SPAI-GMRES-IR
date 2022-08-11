%Example script for comparing GMRES-IR, GMRESIR_SPAI and GMRESIR_NP (with 3
%precisions) 

warning off

%%%%%%%%%% pores_3 %%%%%%%%%%

load('pores_3');
A = Problem.A;

% Find best reordering
[besto, Q] = findbestreordering(A);

fprintf('pores_3, n = %d, nnz(A) = %d, kinf(A) = %0.1e, cond2(A^T) = %0.1e, LU ordering = %s\n', size(A,1), nnz(A), cond(full(A),'inf'), norm(abs(full(A'))*abs(full(inv(A'))),2), besto);


% Run experiments
generateplots(A, Q, 1, 2, 4, 1e-8, .5, 30, 8, 15)
generateplots(A, Q, 1, 2, 4, 1e-8, .4, 30, 8, 15)

fprintf('\n');

%%%%%%%%%% steam1 %%%%%%%%%%

load('steam1');
A = Problem.A;

% Find best reordering
[besto, Q] = findbestreordering(A);

fprintf('steam1, n = %d, nnz(A) = %d, kinf(A) = %0.1e, cond2(A^T) = %0.1e, LU ordering = %s\n', size(A,1), nnz(A), cond(full(A),'inf'), norm(abs(full(A'))*abs(full(inv(A'))),2), besto);


% Run experiments
generateplots(A, Q, 1, 2, 4, 1e-8, .1, 20, 8, 15)
generateplots(A, Q, 1, 2, 4, 1e-8, .2, 20, 8, 15)

fprintf('\n');

%%%%%%%%%% steam3 %%%%%%%%%%

load('steam3');
A = Problem.A;

% Find best reordering
[besto, Q] = findbestreordering(A);

fprintf('steam3, n = %d, nnz(A) = %d, kinf(A) = %0.1e, cond2(A^T) = %0.1e, LU ordering = %s\n', size(A,1), nnz(A), cond(full(A),'inf'), norm(abs(full(A'))*abs(full(inv(A'))),2), besto);


% Run experiments
generateplots(A, Q, 1, 2, 4, 1e-8, .1, 10, 8, 15)
generateplots(A, Q, 1, 2, 4, 1e-8, .5, 10, 8, 15)

fprintf('\n');

%%%%%%%%%% saylr1 %%%%%%%%%%

load('saylr1');
A = Problem.A;

% Find best reordering
[besto, Q] = findbestreordering(A);

fprintf('saylr1, n = %d, nnz(A) = %d, kinf(A) = %0.1e, cond2(A^T) = %0.1e, LU ordering = %s\n', size(A,1), nnz(A), cond(full(A),'inf'), norm(abs(full(A'))*abs(full(inv(A'))),2), besto);


% Run experiments
generateplots(A, Q, 1, 2, 4, 1e-8, .4, 10, 8, 15)
generateplots(A, Q, 1, 2, 4, 1e-8, .3, 30, 8, 15)

fprintf('\n');

%%%%%%%%%% bfwa782 %%%%%%%%%%

load('bfwa782');
A = Problem.A;

% Find best reordering
[besto, Q] = findbestreordering(A);

fprintf('bfwa782, n = %d, nnz(A) = %d, kinf(A) = %0.1e, cond2(A^T) = %0.1e, LU ordering = %s\n', size(A,1), nnz(A), cond(full(A),'inf'), norm(abs(full(A'))*abs(full(inv(A'))),2), besto);

% Run experiments
generateplots(A, Q, 0, 1, 2, 1e-4, .3, 20, 8, 15)
generateplots(A, Q, 0, 1, 2, 1e-4, .5, 20, 8, 15)
generateplots(A, Q, 1, 1, 2, 1e-4, .3, 20, 8, 15)
generateplots(A, Q, 1, 1, 2, 1e-4, .5, 20, 8, 15)

fprintf('\n');

%%%%%%%%%% cage5 %%%%%%%%%%

load('cage5');
A = Problem.A;

% Find best reordering
[besto, Q] = findbestreordering(A);

fprintf('cage5, n = %d, nnz(A) = %d, kinf(A) = %0.1e, cond2(A^T) = %0.1e, LU ordering = %s\n', size(A,1), nnz(A), cond(full(A),'inf'), norm(abs(full(A'))*abs(full(inv(A'))),2), besto);

% Run experiments
generateplots(A, Q, 0, 1, 2, 1e-4, .3, 5, 8, 15)
generateplots(A, Q, 0, 1, 2, 1e-4, .5, 5, 8, 15)
generateplots(A, Q, 1, 1, 2, 1e-4, .3, 5, 8, 15)
generateplots(A, Q, 1, 1, 2, 1e-4, .5, 5, 8, 15)

fprintf('\n');

%%%%%%%%%% gre_115 %%%%%%%%%%

load('gre_115');
A = Problem.A;

% Find best reordering
[besto, Q] = findbestreordering(A);

fprintf('gre_115, n = %d, nnz(A) = %d, kinf(A) = %0.1e, cond2(A^T) = %0.1e, LU ordering = %s\n', size(A,1), nnz(A), cond(full(A),'inf'), norm(abs(full(A'))*abs(full(inv(A'))),2), besto);

% Run experiments
generateplots(A, Q, 0, 1, 2, 1e-4, .3, 15, 8, 15)
generateplots(A, Q, 0, 1, 2, 1e-4, .5, 15, 8, 15)
generateplots(A, Q, 1, 1, 2, 1e-4, .3, 15, 8, 15)
generateplots(A, Q, 1, 1, 2, 1e-4, .5, 15, 8, 15)






