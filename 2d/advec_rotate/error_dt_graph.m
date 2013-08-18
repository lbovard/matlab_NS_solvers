clear all;
clf;

error_inf=dlmread('error_inf');
times = dlmread('times');

loglog(times,error_inf,'.','markersize',20)
hold on
loglog(times,times.^2,'-')
text(1e-3,5e-6,'{\Delta} t^{2}','fontsize',20')
grid on
xlabel {\Delta}t
h=get(gca,'xlabel');
set(h,'FontSize', 20)
ylabel Err_{inf}
h=get(gca,'ylabel');
set(h,'FontSize',20)
title('Convergence 2D Advection w/ leap frog, N=30    ','fontsize',20)