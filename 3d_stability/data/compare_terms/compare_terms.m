function norm_values = compare_terms(kz,F_h,make_plots)

global N k_x k_y k_z k_sq k_inv_sq X Y Fh Re Sc n_k cut u_0 v_0 om_i;

% Constants
N=256; % number of grid point
L=9; %length of box
Fh=F_h;%Froude number
Sc=1; %Schmidt number
Re=10000; %Reynolds number
dx=L/N; %grid step size
dy=dx;
dt=0.0019/3; %time step size
t_final=45;
show_plots=0;	
show_time=0;
makeplots=make_plots;
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

%now get the files
cd ../data

matrixfilename=strcat('k_z.',num2str(kz),'u');
u=dlmread(matrixfilename);
matrixfilename=strcat('k_z.',num2str(kz),'v');
v=dlmread(matrixfilename);
matrixfilename=strcat('k_z.',num2str(kz),'w');
w=dlmread(matrixfilename);
matrixfilename=strcat('k_z.',num2str(kz),'rho');
rho=dlmread(matrixfilename);

%return to previous directory
cd ../compare_terms

%compute useful quantities
u_hat=modfft2(u);
v_hat=modfft2(v);
w_hat=modfft2(w);
rho_hat= modfft2(rho);
omega_x =ifft2(1i*k_y.*w_hat-1i*k_z.*v_hat);
omega_y =ifft2(1i*k_z.*u_hat-1i*k_x.*w_hat);
omega_z =ifft2(1i*k_x.*v_hat-1i*k_y.*u_hat);

%check omega quantities
%norm(omega_z-omega)
%norm(omega)
%norm(omega_z)


%diffusion terms
u_diffusion = real(ifft2(k_sq.*u_hat/Re));
v_diffusion = real(ifft2(k_sq.*v_hat/Re));
w_diffusion = real(ifft2(k_sq.*w_hat/Re));
rho_diffusion=real(ifft2(k_sq.*rho_hat/Re/Sc));
u_diff_norm=norm(u_diffusion);
v_diff_norm=norm(v_diffusion);
w_diff_norm=norm(w_diffusion);
rho_diff_norm=norm(rho_diffusion);

%advection terms
A=v.*om_i+v_0.*omega_y;
B=-u.*om_i-u_0.*omega_z;
C=u_0.*omega_y-v_0.*omega_x;


for i=1:3
	cad{i}=real(ifft2(projection(i,1).*modfft2(v.*om_i)-projection(i,2).*modfft2(u.*om_i)));
	bad{i}=real(ifft2(projection(i,1).*modfft2(v_0.*omega_y)-projection(i,2).*modfft2(u_0.*omega_z)+projection(i,3).*modfft2(C)));
	prho{i}=real(ifft2(projection(i,3).*rho_hat));
end
u_cad=norm(cad{1});
v_cad=norm(cad{2});
w_cad=norm(cad{3});
u_bad=norm(bad{1});
v_bad=norm(bad{2});
w_bad=norm(bad{3});
u_prho_norm=norm(prho{1});
v_prho_norm=norm(prho{2});
w_prho_norm=norm(prho{3});


%other terms
rho_w= real(ifft2(w_hat/Fh^2));
rho_u = real(ifft2(-1i*k_x.*modfft2(u_0.*rho)));
rho_v = real(ifft2(-1i*k_y.*modfft2(v_0.*rho)));
rho_u_norm=norm(rho_u);
rho_v_norm=norm(rho_v);
rho_w_norm=norm(rho_w);


%send back the values
norm_values=[u_diff_norm,v_diff_norm,w_diff_norm,rho_diff_norm,u_cad,v_cad,w_cad,u_bad,v_bad,w_bad,u_prho_norm,v_prho_norm,w_prho_norm,rho_u_norm,rho_v_norm,rho_w_norm];

%now make some plots
if(makeplots==1)
	cd ../plots/term_plots
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'u_diff.ps');
	surf(X,Y,u_diffusion)
	print(h,'-dps',plotfilename);
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'v_diff.ps');
	surf(X,Y,v_diffusion)
	print(h,'-dps',plotfilename);
	clf;
	u=figure;
	plotfilename=strcat('k_z.',num2str(kz),'w_diff.ps');
	surf(X,Y,w_diffusion)
	print(h,'-dps',plotfilename);
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'rho_diff.ps');
	surf(X,Y,rho_diffusion)
	print(h,'-dps',plotfilename);
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'u_cad.ps');
	surf(X,Y,cad{1})
	print(h,'-dps',plotfilename);
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'v_cad.ps');
	surf(X,Y,cad{2})
	print(h,'-dps',plotfilename);
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'w_cad.ps');
	surf(X,Y,cad{3})
	print(h,'-dps',plotfilename);
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'u_bad.ps');
	surf(X,Y,bad{1})
	print(h,'-dps',plotfilename);
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'v_bad.ps');
	surf(X,Y,bad{2})
	print(h,'-dps',plotfilename);
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'w_bad.ps');
	surf(X,Y,bad{3})
	print(h,'-dps',plotfilename);
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'u_prho.ps');
	surf(X,Y,prho{1})
	print(h,'-dps',plotfilename);
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'v_prho.ps');
	surf(X,Y,prho{2})
	print(h,'-dps',plotfilename);
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'w_prho.ps');
	surf(X,Y,prho{3})
	print(h,'-dps',plotfilename);
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'rho_w.ps');
	surf(X,Y,rho_w)
	print(h,'-dps',plotfilename);
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'rho_u.ps');
	surf(X,Y,rho_u)
	print(h,'-dps',plotfilename);
	clf;
	h=figure;
	plotfilename=strcat('k_z.',num2str(kz),'rho_v.ps');
	surf(X,Y,rho_v)
	print(h,'-dps',plotfilename);

	%return to directory
	cd ../../compare_terms
end	
end

