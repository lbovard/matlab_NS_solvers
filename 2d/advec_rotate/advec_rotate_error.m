clear all;

Nvec=2.^(3:8);

for itr=1:length(Nvec)
tic
J=10;
N=Nvec(itr);
h=2*pi/N;
x=linspace(-pi+h,pi,N);
y=x;
[X,Y]=meshgrid(x,y);
dt=h/J;
tot=2*J*N;

%actual functions
c=.5*Y;
d=-.5*X;

%initial wave speeds, approximatetly

c_init=.5;
d_init = 0*0.5;

%init = exp(-(X).^2-(Y).^2);

%define the k matrix
k=repmat([0:N/2-1 0 -N/2+1:-1],N,1)';


u=exp(-10*(X).^2-10*(Y-1).^2);
init=u;
u_old=exp(-10*(X-c_init*dt).^2-10*(Y-d_init*dt-1).^2);

for inner=1:tot
       %check the multiplications so that across rows is x, down columns is y
       fft_u=fft2(u);
       xhat=real(ifft2(1i.*fft_u.*k'));
       yhat=real(ifft2(1i*k.*fft_u));
   
       u_new=u_old - 2*dt*c.*xhat - 2*dt*d.*yhat;
       u_old=u;
       u=u_new;
       %surf(X,Y,u)
       %axis([-pi pi -pi pi -0.1 2])
       %view(2)
       %M(itr)=getframe;
end
error_inf(itr)=norm(u-init,inf);    
toc
end
