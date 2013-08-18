clear all;
clf;
N=2^7;

h=2*pi/N;
%h=2/N;
x=-pi+h*(1:N); %linspace(-pi+h,pi,N);
%x=-1+h*(1:N);
y=x;   
[X,Y]=meshgrid(x,y);

%CLF condition
j=10;
dt=h/j;
tot=2*j*N;
totaltime=0;
nu=h^4;
%nu=10^-4;
end_time=150;
timeitr = floor(end_time/dt);
display=0;

%define the k matrix
k_x=repmat([0:N/2-1 0 -N/2+1:-1],N,1);
k_y=k_x';
k_sq=k_y.^2+k_x.^2;
k_inv=1./k_sq;
%hackish method to find the div by 0, replace
k_zeros=find(k_inv==Inf);
k_inv(k_zeros)=0;

%mean values
use_mean=0;
u_mean=-0.031389066;
v_mean=0;

%co-efficients for the scheme
A=1+nu*dt*k_sq;
B=1-nu*dt*k_sq;
C=B./A;
D=dt*23./(12*A);
E=-16*dt./(12*A);
F=dt*5./(12*A);

%fix the colour scaling
caxis manual;


tic

%omega_old=dlmread('vortex');
%lamb_dipole;
%omega_old=dlmread('dipolev');
%omega_old=exp(-20*(Y).^2)+exp(-20*(X).^2);
epsilon=0.01;	
%omega_old=dlmread('dipole');
%omega_old=exp(-10*(X).^2);
omega_old=2/sqrt(2*pi*epsilon)*exp((-(Y-1.1).^2-(X).^2)/(2*epsilon))-1/sqrt(2*pi*epsilon)*exp((-(Y+1.1).^2-(X).^2)/(2*epsilon));
%omega_old=exp(3*(-(Y).^2-(X-0.5).^2))+exp(3*(-(Y).^2-(X+0.5).^2));
fft_omega_n_2=fft2(omega_old);

%% First Forward Euler
%obtain u,v immediately
if (use_mean==0)
u=real(ifft2(1i*k_y.*fft_omega_n_2.*k_inv));
else
u=1i*k_y.*fft_omega_n_2.*k_inv;
u(1,1)=N^2*u_mean;
u=real(ifft2(u));
end

v=real(ifft2(-1i*k_x.*fft_omega_n_2.*k_inv));

%obtain the omega derivatives
omega_x_n_2=real(ifft2(1i.*k_x.*fft_omega_n_2));
omega_y_n_2=real(ifft2(1i.*k_y.*fft_omega_n_2));
fft_omega_x=dt*fft2(u.*omega_x_n_2);
fft_omega_y=dt*fft2(v.*omega_y_n_2);

fft_omega_n_1=fft_omega_n_2-fft_omega_x-fft_omega_y-nu*dt*k_sq.*fft_omega_n_2;
enstrophy(1)=sum(sum(abs(fft_omega_n_1).^2));
omega_1=real(ifft2(fft_omega_n_1));
energy(1)=sum(sum(u.^2+v.^2));
meanu(1)=mean(u(:));
meanv(1)=mean(v(:));
meanvor(1)=mean(omega_1(:));
[a,b]=max(omega_1(:));
avg_vel(1)=X(b)/dt;

%% Second Forward Euler
%obtain u,v immediately
if (use_mean==0)
u=real(ifft2(1i*k_y.*fft_omega_n_2.*k_inv));
else
u=1i*k_y.*fft_omega_n_2.*k_inv;
u(1,1)=N^2*u_mean;
u=real(ifft2(u));
end
v=real(ifft2(-1i*k_x.*fft_omega_n_1.*k_inv));

%obtain the omega derivatives
omega_x_n_1=real(ifft2(1i*k_x.*fft_omega_n_1));
omega_y_n_1=real(ifft2(1i*k_y.*fft_omega_n_1));
fft_omega_x=dt*fft2(u.*omega_x_n_1);
fft_omega_y=dt*fft2(v.*omega_y_n_1);


fft_omega_n=fft_omega_n_2-fft_omega_x-fft_omega_y-nu*dt*k_sq.*fft_omega_n_1;
enstrophy(2)=sum(sum(abs(fft_omega_n).^2));
omega_2=real(ifft2(fft_omega_n));
energy(2)=sum(sum(u.^2+v.^2));
meanu(2)=mean(u(:));
meanv(2)=mean(v(:));
meanvor(2)=mean(omega_2(:));
[a,b]=max(omega_2(:));
avg_vel(2)=X(b)/(2*dt);

if 1
%% The iteration
for itr=3:timeitr
    %obtain u,v immediately
    if (use_mean==0)
    u=real(ifft2(1i*k_y.*fft_omega_n_2.*k_inv));
    else
    u=1i*k_y.*fft_omega_n_2.*k_inv;
    u(1,1)=N^2*u_mean;
    u=real(ifft2(u));
    end


    v=real(ifft2(-1i*k_x.*fft_omega_n.*k_inv));
    
    %first term
    omega_x_n=real(ifft2(1i*k_x.*fft_omega_n));
    omega_y_n=real(ifft2(1i*k_y.*fft_omega_n));
    a=-fft2(u.*omega_x_n);
    b=-fft2(v.*omega_y_n);
    F1=D.*(a+b);
    
    %second term
    a=-fft2(u.*omega_x_n_1);
    b=-fft2(v.*omega_y_n_1);
    F2=E.*(a+b);
    
    %third term
    a=-fft2(u.*omega_x_n_2);
    b=-fft2(v.*omega_y_n_2);
    F3=F.*(a+b);
    
    
    %time step
    fft_omega_new=C.*fft_omega_n+F1+F2+F3;
    
    %update
    fft_omega_n_2=fft_omega_n_1;
    omega_x_n_2=omega_x_n_1;
    omega_y_n_2=omega_y_n_1;
    
    fft_omega_n_1=fft_omega_n;
    omega_x_n_1=omega_x_n;
    omega_y_n_1=omega_y_n;
    
    fft_omega_n=fft_omega_new;
    
    
    omega=real(ifft2(fft_omega_n));
    enstrophy(itr)=sum(sum(abs(fft_omega_n).^2));
    energy(itr)=sum(sum(u.^2+v.^2));
    meanu(itr)=mean(u(:));
    meanv(itr)=mean(v(:));
    meanvor(itr)=mean(omega(:));
    [a,b]=max(omega(:));
    avg_vel(itr)=X(b)/(itr*dt);

    if display
    if( itr < 1000)
    if mod(itr,25)==0
				surf(X,Y,omega,'EdgeColor','none')
				%axis([-1+h 1 -1+h 1 -1 1])
				axis([-pi+h pi -pi+h pi -15 15])
				view(2);
				M(itr)=getframe;
    end
    else
   if mod(itr,100)==0
				surf(X,Y,omega,'EdgeColor','none')
				axis([-pi+h pi -pi+h pi -15 15])
				view(2);
				M(itr)=getframe;
    end
    end
end
end
end
toc


if 0
 figure(2);
upper_x_limit=length(enstrophy);
upper_y_limit=max(log(enstrophy));
lower_y_limit=min(log(energy));
plot(log(enstrophy))
hold on
plot(log(energy))

axis([0 upper_x_limit lower_y_limit-1 upper_y_limit+1])
figure(3);
surf(X,Y,omega_old,'EdgeColor','none')
axis([-pi+h pi -pi+h pi -10 10])
view(2)
end
