%Example script for comparing GMRES-IR, GMRESIR_SPAI and GMRESIR_NP (with 3 precisions)

load('gre_115');
A = Problem.A;

% Find best reordering
[besto, Q] = findbestreordering(A);
besto

% Run experiments
generateplots(A, Q, 0, 1, 2, 1e-4, .5, 15, 8, 15)

