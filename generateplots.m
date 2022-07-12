%Generate Plot Example script for comparing GMRES-IR, GMRESIR_SPAI and GMRESIR_NP (with 3 precisions)
%   Note: Requires Advanpix multiprecision toolbox
clc;
clear all
data=load('jpwh_991');
A = spconvert(data.Problem.A);
%A=gallery('poisson',10);
n=length(A);
b=ones(n,1);
     %Run GMRES-IR with half,single,double
     fprintf('Running GMRES-IR\n');
     gmresir3(full(A),zeros(n,1),ones(n,1),0,1,2,15,1e-4);
     %Run GMRES-IR with single,single,double
     gmresir3(full(A),zeros(n,1),ones(n,1),1,1,2,15,1e-4);
     %Run GMRES-IR with single,double,quad
     gmresir3(full(A),zeros(n,1),ones(n,1),1,2,4,15,1e-4);
     drawnow
     %Run GMRESIR_SPAI with half,single,double
     fprintf('Running GMRESIR_SPAI\n');
     gmresir3_spai(full(A), zeros(n,1), ones(n,1), .2, 50, 50, 0, 1, 2, 15, 1e-4);
     %Run GMRESIR_SPAI with single,single,double
     gmresir3_spai(full(A), zeros(n,1), ones(n,1), .2, 50, 50, 1, 1, 2, 15, 1e-4);
     %Run GMRESIR_SPAI with single,double,quad
     gmresir3_spai(full(A), zeros(n,1), ones(n,1), .2, 50, 50, 1, 2, 4, 15, 1e-4);
     drawnow
     %Run GMRESIR_NP with half,single,double
     fprintf('Running GMRESIR_np\n');
     gmresir3_np(full(A), zeros(n,1), ones(n,1), 0, 1, 2, 15, 1e-4);
     %Run GMRESIR_NP with single,single,double
     gmresir3_np(full(A), zeros(n,1), ones(n,1), 1, 1, 2, 15, 1e-4);
     %Run GMRESIR_NP with single,double,quad
     gmresir3_np(full(A), zeros(n,1), ones(n,1), 1, 2, 4, 15, 1e-4);
     drawnow
