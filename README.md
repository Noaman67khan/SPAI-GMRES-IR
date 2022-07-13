# SPAI-GMRES-IR
Matlab codes for performing experiments comparing GMRES-IR with LU preconditioning, GMRES-IR with SPAI preconditioning, and GMRES-IR with no preconditioning

## Related publications
Carson, Erin, and Noaman Khan. "Mixed Precision Iterative Refinement with Sparse Approximate Inverse Preconditioning." arXiv preprint arXiv:2202.10204 (2022).

## Included MATLAB files
* **_gmresir3.m_** is a function that performs GMRES-based iterative refinement with LU preconditioning in three precisions, where within GMRES the preconditioner is applied in double the working precision.

* **_sgmresir3.m_** is a function that performs GMRES-based iterative refinement with LU preconditioning in three precisions, where within GMRES a uniform precision is used for all computations (the working precision). 

* **_gmresir3_spai.m_** is a function that performs GMRES-based iterative refinement with SPAI preconditioning in three precisions, where within GMRES the preconditioner is applied in double the working precision.

* **_sgmresir3_spai.m_** is a function that performs GMRES-based iterative refinement with SPAI preconditioning in three precisions, where within GMRES a uniform precision is used for all computations (the working precision). 

* **_gmresir3_np.m_** is a function that performs GMRES-based iterative refinement with no preconditioning in three precisions, where within GMRES the preconditioner is applied in double the working precision.

* **_gmres_hs.m, gmres_sd.m, and gmres_dq.m_** are functions that run left-LU-preconditioned GMRES using precisions half/single, single/double, and double/quad, resp. Application of the preconditioned coefficient matrix to a vector and the preconditioner to the right-hand-side vector are performed in the higher precision; other computations performed all use the lower precision. 

* **_gmres_hh.m, gmres_ss.m, and gmres_dd.m_** are functions that run left-LU-preconditioned GMRES using precisions half, single, and double, resp. All computations are performed in the working precision. 

* **_gmres_spai_XY.m_** are functions that run left-preconditioned GMRES preconditioned using a SPAI preconditioniner using precisions X/Y, as above.

* **_gmres_np_X.m_** are functions that run GMRES with no preconditioning in precision X, as above. 

* **_generateplots.m_** is an example script for experiment on matrix jpwh991, available from SuiteSparse.

## Requirements
* This code requires the Advanpix Multiprecision Computing Toolbox for extended precision computations. 
A free trial of Advanpix is available for download from https://www.advanpix.com/.Advanpix library for quadruple precision. 

* This code also requires functions https://github.com/higham/chop and https://github.com/SrikaraPranesh/Multi_precision_NLA_kernels

