function sigma = simulation(kz,F_h)

global N k_x k_y k_z k_sq k_inv_sq X Y Fh Re Sc n_k cut;

% Constants
N=256; % number of grid point
L=9; %length of box
Fh=F_h;%Froude number
Sc=1; %Schmidt number
Re=16000; %Reynolds number
dx=L/N; %grid step size
dy=dx;
dt=0.0019/2; %time step size
t_final=15;
show_plots=0;	
show_time=0;
num_steps=floor(t_final/dt);
k_z=kz; %vertical wavenumber
T=0:dt:t_final;
n=ceil(2*N/3); %truncation rule
%make truncation rule even
if(mod(n,2)==1)
	n=n-1;
end
%number of terms to be cut off from fft
n_k=(N-n)/2;

% Wave numbers
k_x=2*pi/L*repmat([0:N/2-1 0 -N/2+1:-1],N,1);
k_y=k_x';
k_z=k_z*ones(N,N);
k_sq=k_y.^2+k_x.^2 +k_z.^2;
kzeros=find(k_sq==0);
k_inv_sq=1./k_sq;
k_inv_sq(kzeros)=0;

%truncation
cut = ones(N,N);
cut(:,(N/2-n_k+1):N/2)=0;
cut(:,N/2+2:(N/2+n_k+1))=0;
cut=cut.*cut';

%set stuff up
x=-L/2+dx*(1:N);
y=x';
[X,Y]=meshgrid(x,y);
precalc_projection;
lamb_dipole;
vel_i_hat=initial_velocities;
rho_i_hat=zeros(N,N);


%% Do first step forward euler 
t=0;

P_hat=rho_i_hat.*exp(k_sq*t/Re/Sc);
pr_old=rho_right(vel_i_hat{3},P_hat,t);
P_hat_new=P_hat + dt*pr_old;


for i=1:3
	U_hat{i}=vel_i_hat{i}.*exp(k_sq*t/Re);
end

vr_old=vel_right(U_hat,P_hat,t);
for i=1:3
	U_hat_new{i}=U_hat{i} + dt*vr_old{i};
end

% update scheme
for i=1:3 
	U_hat{i} = U_hat_new{i};
end
P_hat=P_hat_new;

for j=1:num_steps
	t=j*dt;

	pr=rho_right(U_hat{3}.*exp(-k_sq*t/Re),P_hat,t);

	P_hat_new=P_hat+1.5*dt*pr - 0.5*dt*pr_old;

	vr=vel_right(U_hat,P_hat,t);

	for i=1:3
		U_hat_new{i}=U_hat{i} +1.5*dt*vr{i}-0.5*dt*vr_old{i};
	end
	
	%% update scheme 
	P_hat=P_hat_new;
	pr_old=pr;
	display(j)
	for i=1:3
		vr_old{i} =vr{i};
		vel_hat{i}=U_hat{i}.*exp(-k_sq*t/Re);
		U_hat{i} = U_hat_new{i};
	end
	energy(j)=sum(sum(abs(vel_hat{1}).^2+abs(vel_hat{2}).^2+abs(vel_hat{3}).^2));	
	if(show_time==1)
		if(mod(j,100)==1)
			display(t);
		end
	end

        if(show_plots==1)
            if(mod(j,100)==0)
        	surf(X,Y,real(rho),'EdgeColor','none');
        	view(2);
        	M(j)=getframe;
           end
        end
end
 

energy=0.5*log(energy);
m=length(energy);
D=sparse(1:m,1:m,-1/dt,m,m)+sparse(1:m,[2:m 1],1/dt,m,m);
dE=D*energy';	
dE=dE(1:end-1);
sigma=dE(end-1);	
h=figure;
plot(T(1:length(dE)),dE)
xlabel('T')
ylabel('log E')
title_name=strcat('Log Energy vs Time for  k_z=',num2str(kz));
title(title_name);
plotfilename=strcat('k_z.',num2str(kz),'.ps');
print(h,'-dps',plotfilename);
display(plotfilename);
end
