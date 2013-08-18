clear all;
N=3^7;

h=2*pi/N;
x=-pi+h*(1:N); %linspace(-pi+h,pi,N);
y=x;   
[X,Y]=meshgrid(x,y);

%CLF condition
j=1;
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

%define the projection tensor
P_11 = 1-k_x.*k_x.*k_inv;
P_12 = -k_x.*k_y.*k_inv;
P_22 = 1-k_y.*k_y.*k_inv;
P=[P_11,P_12,P_22];

%co-efficients for the scheme
A=1+nu*dt*k_sq;
B=1-nu*dt*k_sq;
c_1=B./A;
c_2=dt*23./(12*A);
c_3=-16*dt./(12*A);
c_4=dt*5./(12*A);


tic

omega_old=dlmread('dipolev');
%omega_old=exp(-20*(Y).^2)+exp(-20*(X).^2);
%%Two Dipoles
%omega_old=exp(-(Y-1).^2-(X).^2)-exp(-(Y+1).^2-(X).^2);
%omega_old=exp(-10*(X).^2);
fft_omega=fft2(omega_old);
%obtain u,v from distribution
u=real(ifft2(-1i*k_y.*fft_omega.*k_inv));
v=real(ifft2(1i*k_x.*fft_omega.*k_inv));
u_hat=fft2(u);
v_hat=fft2(v);

%% First Forward Euler
[F_n_2_hat,G_n_2_hat,omega]=data_projection(u,v,u_hat,v_hat,k_x,k_y,k_inv,k_sq,P_11,P_12,P_22);
if 1
%time-step
u_hat=u_hat-nu*dt*k_sq.*u_hat+dt*F_n_2_hat;
v_hat=v_hat-nu*dt*k_sq.*v_hat+dt*G_n_2_hat;

u=real(ifft2(u_hat));
v=real(ifft2(v_hat));

%compute some quantities
meanu(1)=mean(u(:));
meanv(1)=mean(v(:));
meanvor(1)=mean(omega(:));
energy(1)=sum(sum(u.^2+v.^2));
enstrophy(1)=sum(sum(abs(omega).^2));

%% Secod Forward Euler
[F_n_1_hat,G_n_1_hat,omega]=data_projection(u,v,u_hat,v_hat,k_x,k_y,k_inv,k_sq,P_11,P_12,P_22);

%time-step
u_hat=u_hat-nu*dt*k_sq.*u_hat+dt*F_n_1_hat;
v_hat=v_hat-nu*dt*k_sq.*v_hat+dt*G_n_1_hat;

u=real(ifft2(u_hat));
v=real(ifft2(v_hat));

%compute some quantities
meanu(3)=mean(u(:));
meanv(2)=mean(v(:));
meanvor(2)=mean(omega(:));
energy(2)=sum(sum(u.^2+v.^2));
enstrophy(2)=sum(sum(abs(omega).^2));
%% The iteration
if 1
for itr=2:1000
	%compute the data
	[F_n_hat,G_n_hat,omega]=data_projection(u,v,u_hat,v_hat,k_x,k_y,k_inv,k_sq,P_11,P_12,P_22);

	%time-step
	u_hat=c_1.*u_hat+c_2.*F_n_hat+c_3.*F_n_1_hat+c_4.*F_n_2_hat;
	v_hat=c_1.*v_hat+c_2.*G_n_hat+c_3.*G_n_1_hat+c_4.*G_n_2_hat;

	%update values
	F_n_2_hat=F_n_1_hat;
	G_n_2_hat=G_n_1_hat;
	F_n_1_hat=F_n_hat;
	G_n_1_hat=G_n_hat;

	u=real(ifft2(u_hat));
	v=real(ifft2(v_hat));

	%compute some quantities
	meanu(itr)=mean(u(:));
	meanv(itr)=mean(v(:));
	meanvor(itr)=mean(omega(:));
	energy(itr)=sum(sum(u.^2+v.^2));
	enstrophy(itr)=sum(sum(abs(omega).^2));

	if 1
	%fix the colour scaling
	caxis manual;
	if( itr < 1000)
		if mod(itr,25)==0
			surf(X,Y,omega,'EdgeColor','none')
			axis([-pi+h pi -pi+h pi -3 3])
			caxis([-1 1])
			view(2);
			M(itr)=getframe;
		end
		else
		if mod(itr,100)==0
			surf(X,Y,omega,'EdgeColor','none')
			axis([-pi+h pi -pi+h pi -3 3])
			caxis([-1 1])
			view(2);
			M(itr)=getframe;
		end
	end
	end
end
end
toc
end
