espai=[0.01,0.1,0.2,0.3,0.4,0.5];
alpha=50;
beta=50;
kinfAts = [];
num = 1;
for i = 1:numel(espai)
    A = gallery('poisson',16);
    %data=load('saylr1.mtx');
    %A = spconvert(data);
    M = spai_ss(A',espai(i),alpha,beta);
    M=M';
    Ats = double(mp(double(M),34))*mp(double(A),34);
    kinfAts1 = cond(mp(double(Ats),34),'inf');
    kinfAts =[kinfAts kinfAts1];
end
kinfAtd = [];
num = 1;
for i = 1:numel(espai)
    A = gallery('poisson',16);
    %data=load('saylr1.mtx');
    %A = spconvert(data);
    M = spai_dd(A',espai(i),alpha,beta);
    M=M';
    Atd = double(mp(double(M),34))*mp(double(A),34);
    kinfAtd1 = cond(mp(double(Atd),34),'inf');
    kinfAtd =[kinfAtd kinfAtd1];
end
fig1 = figure();
semilogy(espai,kinfAtd, '-rx');
hold on
semilogy(espai,kinfAts, '-bo');
set(gca,'FontSize',8)
xlabel ('$\varepsilon$','Interpreter','latex')
ylabel ('$\kappa_{\infty}(\tilde{A})$','Interpreter','latex')
title('Saylr1', 'Interpreter','latex')
