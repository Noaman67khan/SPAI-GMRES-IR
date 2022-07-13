function x = gmresir3_np(A, b, precf, precw, precr, iter_max, gtol)
%GMRESIR3_NP  GMRES-based iterative refinement in three precisions with no
%             preconditoning. 
%     Solves Ax = b using gmres-based
%     iterative refinement with at most iter_max ref. steps and GMRES convergence
%     tolerance gtol, with
%     M computed in precision precf:
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
%
% Note: requires Advanpix multiprecision toolbox, 
% chop library (https://github.com/higham/chop), and 
% https://github.com/SrikaraPranesh/Multi_precision_NLA_kernels

if precf ~=0 && precf ~=1 && precf ~= 2, error('precf should be 0, 1 or 2'), end
if precw ~=0 && precw ~=1 && precw ~= 2, error('precw should be 0, 1 or 2'), end
if precr ~=1 && precr ~= 2 && precr ~= 4, error('precr should be 1, 2, or 4'), end

n = length(A);

if precf == 1
    %fprintf('**** Factorization precision is single.\n')
    ufs = 'S';
elseif precf == 2
    %fprintf('**** Factorization precision is double.\n')
    ufs = 'D';
else
    %fprintf('**** Factorization precision is half.\n')
    ufs = 'H';
end

if precw == 0
    %fprintf('**** Working precision is half.\n')
    uws = 'H';
    A = chop(A);
    b = chop(b);
    u = eps(chop(1));
elseif precw == 2
    %fprintf('**** Working precision is double.\n')
    uws = 'D';
    A = double(A);
    b = double(b);
    u = eps('double');
else
    %fprintf('**** Working precision is single.\n')
    uws = 'S';
    A = single(A);
    b = single(b);
    u = eps('single');
end

if precr == 1
    %fprintf('**** Residual precision is single.\n')
    urs = 'S';
elseif precr == 2
    %fprintf('**** Residual precision is double.\n')
    urs = 'D';
else
    %fprintf('**** Residual precision is quad.\n')
    urs = 'Q';
    mp.Digits(34);
end

xact = double(mp(double(A),34)\mp(double(b),34));  


if precf == 1
    x = single(b);
elseif precf == 2
    x = double(b);
else
    x = chop(b);
end

%Compute condition number of A, of preconditioned system At, cond(A), and
%cond(A,x) for the exact solution
kinfA = cond(mp(double(A),34),'inf');
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
         [d, err, its, ~] = gmres_np_h( A, chop(zeros(n,1)), chop(rd1),n, 1, gtol);
    elseif precw == 2
        [d, err, its, ~] = gmres_np_d( A, zeros(n,1), double(rd1), n, 1, gtol);
    else
          [d, err, its, ~] = gmres_np_s( A, single(zeros(n,1)), single(rd1), n, 1, gtol);
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


%Generate plots

%Create ferr, nbe, cbe plot
fig1 = figure();
semilogy(0:iter-1, ferr, '-rx');
hold on
semilogy(0:iter-1, nbe, '-bo');
hold on
semilogy(0:iter-1, cbe, '-gv');
hold on
semilogy(0:iter-1, double(u)*ones(iter,1), '--k');

if (nargin==10)
xlim([0 lim_num])
ylim([10.^(-30) 10.^(1)]);
xx = lim_num-numel(ferr)+2;
hold on
semilogy(numel(nbe)-1:lim_num, double(u)*ones(xx,1), '--k');
end
%Ensure only integers labeled on x axis
atm = get(gca,'xticklabels');
m = str2double(atm);
xlab = [];
num = 1;
for i = 1:numel(m)
    if ceil(m(i)) == m(i)
        xlab(num) = m(i);
        num = num + 1;
    end
end
set(gca,'xticklabels',xlab);
set(gca,'xtick',xlab);
xlabel({'refinement step'},'Interpreter','latex');

str_e = sprintf('%0.1e',kinfA);
%str_eps = sprintf('%0.1f',espai);
%iter = sprintf('GMRES its = %s\n', num2str(gmresits));
tt = strcat('GMRES-IR Without P,  $$\, \kappa_{\infty}(A) = $$ ',str_e,'');    
title(tt,'Interpreter','latex');

set(gca,'FontSize',15)
a = get(gca,'Children');
set(a,'LineWidth',1);
set(a,'MarkerSize',10);

h = legend('ferr','nbe','cbe');
set(h,'Interpreter','latex');

% %Create phi plot
%fig2 = figure();
%semilogy(0:iter-2, lim, '-cx');
%hold on
%semilogy(0:iter-2, lim2, '-+','Color',[1 0.600000023841858 0.200000002980232]);
%hold on
%semilogy(0:iter-2, etai, '-mo');
%hold on
%semilogy(0:iter-2, phi, '-kv');
%hold on
%semilogy(0:iter-1, ones(iter,1), '--k');

%Use same x labels as error plot
set(gca,'xticklabels',xlab);
set(gca,'xtick',xlab);
xlabel({'refinement step'},'Interpreter','latex');

title(tt,'Interpreter','latex');

%h = legend('$u_s \Vert E_i \Vert_\infty$','$\phi_i$');
%set(h,'Interpreter','latex');
% if ~isempty(savename)
%     saveas(gcf, strcat(savename,'.pdf'));
% end
% if ~isempty(savename)
%     saveas(gcf, strcat(savename,'.pdf'));
% end
[L,U]=lu(A);
fprintf('nnz(A) = %d, nnz(L+U) = %d, nnz(inv(A)) = %d\n', nnz(A), nnz(L+U), nnz(inv(A)));
fprintf('GMRES its = %s\n', num2str(gmresits));

end
