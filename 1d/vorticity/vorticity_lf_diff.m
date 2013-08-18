clear all;
clf;
N=2^6;

h=2*pi/N;
x=-pi+h*(1:N); %linspace(-pi+h,pi,N);

%CLF condition
j=100 ;
dt=h/j;
tot=2*j*N;
totaltime=0;
nu=h^2;

%define the k matrix
k_x=[0:N/2-1 0 -N/2+1:-1];
k_sq=k_x.^2;
k_inv=1./k_sq;
%hackish method to find the div by 0, replace
k_zeros=find(k_inv==Inf);
k_inv(k_zeros)=0;


tic

%omega_old=exp(-10*(Y).^2);
omega_old=exp(-10*(x).^2);
%omega_old=sin(X).*sin(Y);
%omega_old=exp(-(Y).^2-(X-1).^2)-exp(-(Y).^2-(X+1).^2);
fft_omega=fft(omega_old);
%obtain u,v immediately
%u=ifft2(1i*k_y.*fft_omega.*k_inv);
u=1;

enstrophy(1)=sum(sum(abs(fft_omega).^2));

%omega derivatives in x and y
xhat=real(ifft(1i.*k_x.*fft_omega));

%omega laplacian 
lapomega=nu*real(ifft(-k_sq.*fft_omega));

omega=omega_old - dt*(u.*xhat)+dt*lapomega;
if 1
for itr=2:5000
    fft_omega=fft(omega);
    xhat=real(ifft(1i.*fft_omega.*k_x));
    lapomega=nu*real(ifft(-k_sq.*fft_omega));
    omega_new=omega_old - 2*dt*u.*xhat+2*dt*lapomega;
    omega_old=omega;
    omega=omega_new;
    enstrophy(itr)=sum(sum(abs(fft_omega).^2));
    if 1
    if mod(itr,100)==0
     plot(x,omega);
     axis([-pi+h pi -2 2])
     M(itr)=getframe;
    end
    end
end
end
toc
upper_x_limit=length(enstrophy);
upper_y_limit=max(log(enstrophy));
plot(log(enstrophy))
axis([0 upper_x_limit upper_y_limit-1 upper_y_limit+1])