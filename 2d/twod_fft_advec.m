clear all;

N=50;
h=2*pi/N;
x=linspace(-pi,pi,N);
y=x;
[X,Y]=meshgrid(x,y);
dt=h/8;
tot=2*pi/dt;
%actual functions
c=.5*Y;
d=-.5*X;

%initial wave speeds, approximatetly
c_init=.5;
d_init = 0.5;

%init = exp(-(X).^2-(Y).^2);

%define the k matrix
k=repmat([0:N/2-1 0 -N/2+1:-1],N,1)';


tic

u=exp(-10*(X-pi/2).^2-10*(Y-pi/2).^2);
u_old=exp(-10*(X-c_init*dt-pi/2).^2-10*(Y-d_init*dt-pi/2).^2);

for itr=1:2*tot
       %check the multiplications so that across rows is x, down columns is y
       xhat=real(ifft2(1i*fft2(u).*k'));
       yhat=real(ifft2(1i*k.*fft2(u)));
   
       u_new=u_old - 2*dt*c.*xhat - 2*dt*d.*yhat;
       u_old=u;
       u=u_new;
       surf(X,Y,u)
       axis([-pi pi -pi pi -0.1 2])
       view(2)
       M(itr)=getframe;
end
%norm(u-init,inf);
toc
