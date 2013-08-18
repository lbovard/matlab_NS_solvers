clear all;
N=2^7;

h=2*pi/N;
x=-pi+h*(1:N); %linspace(-pi+h,pi,N);
y=x;   
[X,Y]=meshgrid(x,y);

%CLF condition
j=100;
dt=h/j;
tot=2*j*N;
totaltime=0;
nu=h^2;

%define the k matrix
k_x=repmat([0:N/2-1 0 -N/2+1:-1],N,1);
k_y=k_x';
k_sq=k_y.^2+k_x.^2;
k_inv=1./k_sq;
%hackish method to find the div by 0, replace
k_zeros=find(k_inv==Inf);
k_inv(k_zeros)=0;

%co-efficients for the scheme
A=1+nu*dt*k_sq;
B=1-nu*dt*k_sq;
BdivA=A./B;
Binv=1./B;

tic

%omega_old=exp(-10*(Y).^2);
%omega_old=exp(-10*(Y).^2)+exp(-10*(X).^2);
%omega_old=sin(X).*sin(Y);
epsilon=0.1;
%omega_old=dlmread('dipole');
%omega_old=exp(-10*(X).^2);
omega_old=1/sqrt(2*pi*epsilon)*exp((-(Y-0.5).^2-(X).^2)/(2*epsilon))+1/sqrt(2*pi*epsilon)*exp((-(Y+0.5).^2-(X).^2)/(2*epsilon));
%omega_old=exp(-(Y).^2-(X-1).^2)+exp(-(Y).^2-(X+1).^2);
fft_omega_old=fft2(omega_old);
%obtain u,v immediately
u=real(ifft2(1i*k_y.*fft_omega_old.*k_inv));
v=real(ifft2(-1i*k_x.*fft_omega_old.*k_inv));

enstrophy(1)=sum(sum(abs(fft_omega_old).^2));

C=fft2(u.*real(ifft2(1i*k_x.*fft_omega_old)));
D=fft2(v.*real(ifft2(1i*k_y.*fft_omega_old)));

fft_omega=fft_omega_old-dt*C-dt*D;

omega=real(ifft2(fft_omega));
if 1
for itr=2:1500
    %obtain u,v immediately
    u=real(ifft2(1i*k_y.*fft_omega.*k_inv));
    v=real(ifft2(-1i*k_x.*fft_omega.*k_inv));

    C=fft2(u.*real(ifft2(1i*k_x.*fft_omega)));
    D=fft2(v.*real(ifft2(1i*k_y.*fft_omega)));
    
    fft_omega_new=BdivA.*fft_omega_old-2*dt*Binv.*C-2*dt*Binv.*D;
    fft_omega_old=fft_omega;
    fft_omega=fft_omega_new;
    
    omega=real(ifft2(fft_omega));
    enstrophy(itr)=sum(sum(abs(fft_omega).^2));
    if 1
    if mod(itr,100)==0
     surf(X,Y,omega,'EdgeColor','none')
     axis([-pi+h pi -pi+h pi -2 10])
     view(2)
     M(itr)=getframe;
    end
    end
end
end
toc
upper_x_limit=length(enstrophy);
upper_y_limit=max(log(enstrophy));
lower_y_limit=min(log(enstrophy));
plot(log(enstrophy))

axis([0 upper_x_limit lower_y_limit-1 upper_y_limit+1])

