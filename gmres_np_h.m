function [x, error, its, flag] = gmres_np_h( A, x, b, restrt, max_it, tol)
%GMRES_NP_H   Unpreconditioned GMRES in half precision
%   Solves Ax=b 
%   using the Generalized Minimal residual ( GMRES ) method.
%   Currently uses relative residual norm to check for convergence 
%   (same as Matlab GMRES)
%   Half precision used throughout.
%
%   input   A        REAL nonsymmetric positive definite matrix
%           x        REAL initial guess vector
%           b        REAL right hand side vector
%           restrt   INTEGER number of iterations between restarts
%           max_it   INTEGER maximum number of iterations
%           tol      REAL error tolerance
%
%   output  x        REAL solution vector
%           error    REAL error norm
%           its     INTEGER number of (inner) iterations performed
%           flag     INTEGER: 0 = solution found to tolerance
%                             1 = no convergence given max_it
%


flag = 0;
its = 0;

%Ensure half working precision
A = chop(A);
b = chop(b);
x = chop(x);



rtmp = chop(chop(b)-chop(chop(A)*chop(x)));
r = chop(rtmp);

bnrm2 = norm(r );
if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end
error=[];

error(1) = norm( r ) / bnrm2;
if ( error(1) < tol ) return, end

[n,~] = size(A);                                  % initialize workspace
m = restrt;
V = chop(zeros(n,m+1));
H = chop(zeros(m+1,m));
cs = chop(zeros(m,1));
sn = chop(zeros(m,1));
e1    = chop(zeros(n,1));
e1(1) = chop(1.0);

for iter = 1:max_it,                              % begin iteration
    rtmp = chop(chop(b)-chop(chop(A)*chop(x)));
    r = chop(rtmp);
    
    V(:,1) = chop(chop(r) / chop(norm( chop(r) )));
    s = chop(norm( chop(r) )*chop(e1));
    for i = 1:m,                     % construct orthonormal basis via GS
        its = its+1;
        vcur = chop(V(:,i));      
        
        vcur = chop(A)*chop(vcur);
        
        w = chop(vcur);
      
        for k = 1:i,
            H(k,i)= chop(chop(w')*chop(V(:,k)));
            w = chop(chop(w) - chop(chop(H(k,i))*chop(V(:,k))));
        end
        H(i+1,i) = norm( w );
        V(:,i+1) = chop(chop(w) / chop(H(i+1,i)));
        for k = 1:i-1                              % apply Givens rotation
            temp     =  chop(chop(chop(cs(k))*chop(H(k,i))) + chop(chop(sn(k))*chop(H(k+1,i))));
            H(k+1,i) = chop(chop(chop(-sn(k))*chop(H(k,i))) + chop(chop(cs(k))*chop(H(k+1,i))));
            H(k,i)   = temp;
        end
        [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); % form i-th rotation matrix
        temp   = chop(chop(cs(i))*chop(s(i)));                        % approximate residual norm
        s(i+1) = chop(chop(-sn(i))*chop(s(i)));
        s(i)   = temp;
        H(i,i) = chop(chop(chop(cs(i))*chop(H(i,i))) + chop(chop(sn(i))*chop(H(i+1,i))));
        H(i+1,i) = 0.0;
        error((iter-1)*m+i+1)  = chop(chop(abs(s(i+1))) / chop(bnrm2));
        if ( error((iter-1)*m+i+1) <= tol ),                        % update approximation
            y = chop(chop(H(1:i,1:i)) \ chop(s(1:i)));                 % and exit
            addvec = chop(chop(V(:,1:i))*chop(y));
            x = chop(chop(x) + chop(addvec));
            break;
        end
    end
    
    if ( error(end) <= tol ), break, end
    y = chop(chop(H(1:m,1:m)) \ chop(s(1:m)));
    addvec = chop(chop(V(:,1:m))*chop(y));
    x = chop(chop(x) + chop(addvec));                            % update approximation
    rtmp = chop(chop(b)-chop(chop(A)*chop(x)));
    r = chop(r);
    s(i+1) = chop(norm(r));
    error = [error, chop(chop(s(i+1)) / chop(bnrm2))];                        % check convergence
    if ( error(end) <= tol ), break, end;
end

if ( error(end) > tol ) flag = 1; end;                 % converged




function [ c, s ] = rotmat( a, b )
%
% Compute the Givens rotation matrix parameters for a and b.
%
if ( b == 0.0 ),
    c = 1.0;
    s = 0.0;
elseif ( abs(b) > abs(a) ),
    temp = a / b;
    temp2 = temp*temp; 
    s = 1.0 / chop(sqrt(chop( 1.0 + temp2 )));
    c = temp * s;
else
    temp = b / a;
    temp2 = temp*temp;
    c = 1.0 / chop(sqrt(chop( 1.0 + temp2 )));
    s = temp * c;
end
