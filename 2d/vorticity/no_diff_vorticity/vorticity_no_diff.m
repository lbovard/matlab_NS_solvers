clear all;
N=2^7;

h=2*pi/N;
alpha=0.15;
x=-pi+h*(1:N); %linspace(-pi+h,pi,N);
y=x;   
[X,Y]=meshgrid(x,y);
u_mean=0.9204;
v_mean=-3.0921e-06;

%CLF condition
j=10;
dt=h/j;
tot=2*j*N;
totaltime=0;
nu=0;
%lamb_dipole;
%define the k matrix
k_x=repmat([0:N/2-1 0 -N/2+1:-1],N,1);
k_y=k_x';
k_sq=k_y.^2+k_x.^2;
k_inv=1./k_sq;
%hackish method to find the div by 0, replace
k_zeros=find(k_inv==Inf);
k_inv(k_zeros)=0;


tic

omega_old=dlmread('dipolev');
%omega_old=exp(-10*((Y).^2+(X-1).^2))+exp(-10*((Y).^2+(X+1).^2));
fft_omega=fft2(omega_old);
%obtain u,v immediately
u=1i*k_y.*fftnZZ_omega.*k_inv;
u(1,1)=N^2*u_mean;
u=ifft2(u);
v=ifft2(-1i*k_x.*fft_omega.*k_inv);
%u=1;v=1;

enstrophy(1)=sum(sum(abs(fft_omega).^2));

%omega derivatives in x and y
xhat=real(ifft2(1i.*k_x.*fft_omega));
yhat=real(ifft2(1i*k_y.*fft_omega));

%omega laplacian 
lapomega=nu*real(ifft2(-k_sq.*fft_omega));

omega=omega_old - dt*(u.*xhat +v.*yhat)+dt*lapomega;
if 1
for itr=2:1000
    fft_omega=fft2(omega);
    %obtain u,v immediately
    %u=ifft2(1i*k_y.*fft_omega.*k_inv);
    u=1i*k_y.*fft_omega.*k_inv;
    u(1,1)=N^2*u_mean;
    u=ifft2(u);
    v=ifft2(-1i*k_x.*fft_omega.*k_inv);
    xhat=real(ifft2(1i.*fft_omega.*k_x));
    yhat=real(ifft2(1i*k_y.*fft_omega));   
    lapomega=nu*real(ifft2(-k_sq.*fft_omega));
    omega_new=omega_old - 2*dt*u.*xhat - 2*dt*v.*yhat+2*dt*lapomega;
    omega_old=omega;
    omega=omega_new;
    enstrophy(itr)=sum(sum(abs(fft_omega).^2));
    if 1
    if mod(itr,1)==0
     surf(X,Y,omega,'EdgeColor','none')
     axis([-pi+h pi -pi+h pi -20 20])
     view(2)
     M(itr)=getframe;
    end
    end
end
end
toc
upper_x_limit=length(enstrophy);
upper_y_limit=max(log(enstrophy));
%plot(log(enstrophy))
%axis([0 upper_x_limit upper_y_limit-1 upper_y_limit+1])

