function x = sgmresir3(A,Q,b,precf,precw,precr,iter_max, gtol)
%SGMRESIR   GMRES-based iterative refinement in three precisions using LU preconditioning.
%     Solves Ax = b using GMRES-based
%     iterative refinement (with at most iter_max ref. steps), using double 
%     the working precision for applying the preconditioned matrix, with
%     LU factors computed in precision precf:
%       * half if precf = 0,
%       * single if precf = 1,
%       * double if precf = 2,
%     working precision precw:
%       * half if precw = 0,
%       * single if precw = 1,
%       * double if precw = 2,
%     and residuals computed at precision precr:
%       * single if precr = 1,
%       * double if precr = 2,
%       * quad if precr = 4
%     The working precision is used throughout the calls to GMRES
%
% Note: requires Advanpix multiprecision toolbox, 
% chop library (https://github.com/higham/chop), and 
% https://github.com/SrikaraPranesh/Multi_precision_NLA_kernels

if precf ~=0 && precf ~=1 && precf ~= 2, error('precf should be 0, 1 or 2'), end
if precw ~=0 && precw ~=1 && precw ~= 2, error('precw should be 0, 1 or 2'), end
if precr ~=1 && precr ~= 2 && precr ~= 4, error('precr should be 1, 2, or 4'), end

n = length(A);

if precf == 1
   % fprintf('**** Factorization precision is single.\n')
    ufs = 'single';
elseif precf == 2
   % fprintf('**** Factorization precision is double.\n')
    ufs = 'double';
else
   % fprintf('**** Factorization precision is half.\n')
    ufs = 'half';
    fp.format = 'h';
    chop([],fp);
    
end

if precw == 0
   % fprintf('**** Working precision is half.\n')
    uws = 'half';
    uws1 = 'h';
    u = float_params('h');
elseif precw == 2
   % fprintf('**** Working precision is double.\n')
    uws = 'double';
    A = double(A);
    b = double(b);
    u = eps('double');
else
  %  fprintf('**** Working precision is single.\n')
    uws = 'single';
    A = single(A);
    b = single(b);
    u = eps('single');
end

if precr == 1
   % fprintf('**** Residual precision is single.\n')
    urs = 'single';
elseif precr == 2
   % fprintf('**** Residual precision is double.\n')
    urs = 'double';
else
  %  fprintf('**** Residual precision is quad.\n')
    urs = 'quad';
    mp.Digits(34);
end

xact = double(mp(double(A),34)\mp(double(b),34));


[uh,xmins,xmin,xmax] = float_params(ufs);



%Compute LU factorization
if precf == 1
    [L,U,P] = lu(single(Q*A*Q'));
    nnzLU = nnz(L+U);
    J = eye(n);
    LL = single(double(Q')*double(P')*double(L));
    U = single(U*Q);  
    x =  U\(LL\(single(b)) );
    %x =  U\(L\(P*single(b)) );

elseif precf == 2
    [L,U,P] = lu(double(A));
    nnzLU = nnz(L+U);
    LL = double(double(Q')*double(P')*double(L));
    U = U*Q;
    x =  U\(LL\(double(b)) );
    %x =  U\(L\(P*double(b)) );
else
    Ah = chop(A);
    [L,U,pp] = lutx_chop(Q*Ah*Q');
    nnzLU = nnz(L+U);
    I = chop(eye(n)); P = I(pp,:);
    LL = chop(Q'*P'*L);
    U = chop(U);
    
    % If L or U factors contain NAN or INF, try again with scaling
    if ( sum(sum(isinf(single(LL))))>0 || sum(sum(isnan(single(LL))))>0 || sum(sum(isinf(single(U))))>0 || sum(sum(isnan(single(U))))>0 )
        
        [Ah,R,C] = scale_diag_2side(Q*A*Q');
        mu = (0.1)*xmax;
        Ah = mu*Ah;
        
        Ah = chop(Ah);
        [L,U,p] = lutx_chop(Ah);
        
        I = chop(eye(n)); P = I(p,:);
        LL = Q'*P'*L; 
        LL = (1/mu)*diag(1./diag(R))*(LL);
        U = (U*Q)*diag(1./diag(C));
        x = U\(LL\b);
    else
        t1 = lp_matvec(Q,chop(b));
        t1 = lp_matvec(P,chop(t1));
        t1 = trisol(L,b);
        y = trisol(U,t1);
        x = lp_matvec(Q', chop(y));
        U = chop(U*Q);
    end
    
end

%Compute condition number of A, of preconditioned system At, cond(A), and
%cond(A,x) for the exact solution
kinfA = cond(mp(double(A),34),'inf');
At = double(mp(double(U),34)\(mp(double(LL),34)\( mp(double(A),34))));
kinfAt = cond(mp(double(At),34),'inf');
condAx = norm(abs(inv(mp(double(A),34)))*abs(mp(double(A),34))*abs(xact),inf)/norm(xact,inf);
condA = norm(abs(inv(mp(double(A),34)))*abs(mp(double(A),34)),'inf');

%Note: when kinf(A) is large, the initial solution x can have 'Inf's in it
%If so, default to using 0 as initial solution
if sum(isinf(single(x)))>0
    x =  zeros(size(b,1),1);
    fprintf('**** Warning: x0 contains Inf. Using 0 vector as initial solution.\n')
end

%Store initial solution in working precision
if precw == 0
    x = chop(x);
elseif precw == 2
    x = double(x);
else
    x = single(x);
end

cged = false;
iter = 0; dx = 0; rd = 0;

%Array to store total number of gmres iterations in each ref step
gmresits = [];

%Array to store final relative (preconditioned) residual norm in gmres
gmreserr = [];

while ~cged
    
   %Compute size of errors, quantities in bounds
    ferr(iter+1) = double(norm(mp(double(x),34)-mp(xact,34),'inf')/norm(mp(xact,34),'inf'));
    mu(iter+1) = norm(double(A)*(mp(double(x),34)-mp(xact,34)),'inf')/(norm(mp(double(A),34),'inf')*norm(mp(double(x),34)-mp(xact,34),'inf')); 
    res = double(b) - double(A)*double(x);
    nbe(iter+1) = double(norm(mp(res,34),'inf')/(norm(mp(double(A),34),'inf')*norm(mp(double(x),34),'inf')+ norm(mp(double(b),34),'inf')));
    temp = double( abs(mp(res,34)) ./ (abs(mp(double(A),34))*abs(mp(double(x),34)) + abs(mp(double(b),34))) );
    temp(isnan(temp)) = 0; % Set 0/0 to 0.
    cbe(iter+1) = max(temp);
    
    iter = iter + 1;
    if iter > iter_max, break, end
    
    %Check convergence
    if max([ferr(iter) nbe(iter) cbe(iter)]) <= u, break, end
    
    %Compute residual vector
    if precr == 1
        rd = single(b) - single(A)*single(x);
    elseif precr == 2
        rd = double(b) - double(A)*double(x);
    else
        rd = mp(double(b),34) - mp(double(A),34)*mp(double(x),34);
    end
    
    %Scale residual vector
    norm_rd = norm(rd,inf);
    rd1 = rd/norm_rd;
    
    %Call GMRES to solve for correction term
    if precw == 0
        [d, err, its, ~] = gmres_hh( A, chop(zeros(n,1)), chop(rd1), LL, U, n, 1, gtol);
    elseif precw == 2
        [d, err, its, ~] = gmres_dd( A, zeros(n,1), double(rd1), LL, U, n, 1, gtol);
    else
        [d, err, its, ~] = gmres_ss( A, single(zeros(n,1)), single(rd1), LL, U, n, 1, gtol);
    end
    
    
    %Compute quantities in bounds for plotting
    lim(iter) = double( 2*u*cond(mp(double(A),34),'inf')*mu(iter));
    lim2(iter) = double(2*u*condA);  
    dact = mp(double(A),34)\mp(double(rd1),34);
    etai(iter) = norm(double(mp(double(d),34)-dact),'inf')/norm(dact,'inf');   
    phi(iter) = min(lim(iter),lim2(iter))+etai(iter);
    
    %Record number of iterations gmres took
    gmresits = [gmresits,its];
    
    %Record final relative (preconditioned) residual norm in GMRES
    gmreserr = [gmreserr,err(end)];
    
    %Record relative (preconditioned) residual norm in each iteration of
    %GMRES (so we can look at convergence trajectories if need be)
    gmreserrvec{iter} = err;
    
    xold = x;
    
    %Update solution
    if precw == 0
        x = x + chop(norm_rd)*chop(d);
    elseif precw == 2
        x = x + norm_rd*double(d);
    else
        x = x + single(norm_rd)*single(d);
    end
    dx = norm(x-xold,'inf')/norm(x,'inf');
    
    %Check if dx contains infs, nans, or is 0
    if dx == Inf || isnan(double(dx))
        plt = 0;
        break;
    end
    
end

% %Generate plots
% 
% %Create ferr, nbe, cbe plot
% fig1 = figure();
% semilogy(0:iter-1, ferr, '-rx');
% hold on
% semilogy(0:iter-1, nbe, '-bo');
% hold on
% semilogy(0:iter-1, cbe, '-gv');
% hold on
% semilogy(0:iter-1, double(u)*ones(iter,1), '--k');
% 
% if (nargin==10)
% xlim([0 lim_num])
% ylim([10.^(-30) 10.^(1)]);
% xx = lim_num-numel(ferr)+2;
% hold on
% semilogy(numel(nbe)-1:lim_num, double(u)*ones(xx,1), '--k');
% end
% %Ensure only integers labeled on x axis
% atm = get(gca,'xticklabels');
% m = str2double(atm);
% xlab = [];
% num = 1;
% for i = 1:numel(m)
%     if ceil(m(i)) == m(i)
%         xlab(num) = m(i);
%         num = num + 1;
%     end
% end
% set(gca,'xticklabels',xlab);
% set(gca,'xtick',xlab);
% xlabel({'refinement step'},'Interpreter','latex');
% 
% str_e = sprintf('%0.1e',kinfA);
 str_a = sprintf('%0.1e',kinfAt);
% %str_eps = sprintf('%0.1f',espai);
% %iter = sprintf('GMRES its = %s\n', num2str(gmresits));
% tt = strcat('GMRES-IR,  $$\, \kappa_{\infty}(A) = $$ ',str_e,'');    
% title(tt,'Interpreter','latex');
% 
% set(gca,'FontSize',15)
% a = get(gca,'Children');
% set(a,'LineWidth',1);
% set(a,'MarkerSize',10);
% 
% h = legend('ferr','nbe','cbe');
% set(h,'Interpreter','latex');
% 
% % %Create phi plot
% %fig2 = figure();
% %semilogy(0:iter-2, lim, '-cx');
% %hold on
% %semilogy(0:iter-2, lim2, '-+','Color',[1 0.600000023841858 0.200000002980232]);
% %hold on
% %semilogy(0:iter-2, etai, '-mo');
% %hold on
% %semilogy(0:iter-2, phi, '-kv');
% %hold on
% %semilogy(0:iter-1, ones(iter,1), '--k');
% 
% %Use same x labels as error plot
% set(gca,'xticklabels',xlab);
% set(gca,'xtick',xlab);
% xlabel({'refinement step'},'Interpreter','latex');
% 
% title(tt,'Interpreter','latex');
% 
% %h = legend('$u_s \Vert E_i \Vert_\infty$','$\phi_i$');
% %set(h,'Interpreter','latex');
% % if ~isempty(savename)
% %     saveas(gcf, strcat(savename,'.pdf'));
% % end
% % if ~isempty(savename)
% %     saveas(gcf, strcat(savename,'.pdf'));
% % end

fprintf('nnz(A) = %d, nnz(L+U) = %d, nnz(inv(A)) = %d, kinf(At) = %s\n', nnz(A), nnzLU, nnz(inv(A)), str_a);
fprintf('GMRES its = %s(%s)\n', num2str(sum(gmresits)),num2str(gmresits));

end
