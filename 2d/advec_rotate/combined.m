clear all;

Nvec=2.^(4);
for outer=1:length(Nvec)
N=Nvec(outer);
h=2*pi/N;
x=-pi+h*[1:N]; %linspace(-pi+h,pi,N);
y=x;
[X,Y]=meshgrid(x,y);
%used for CLF condition
j=10;
dt=h/j;
tot=2*j*N;
totaltime=0;
%window=(tanh((X+pi-h)/0.1)-tanh((X-pi)/0.1)-1).*(tanh((Y+pi-h)/0.1)-tanh((Y-pi)/0.1)-1);
%advection co-efficients 
c=0.5*ones(N);
d=zeros(N);
c_one=0.5;

%define the k matrix
k=repmat([0:N/2-1 0 -N/2+1:-1],N,1)';
k_one=[0:N/2-1 0 -N/2+1:-1];

tic
u=exp(-10*(X).^2);
u_one=exp(-10*(x).^2);
%u=sin(X).*cos(Y);
init=u;
init_one=u_one;

for itr=1:tot
       %first RK4 term q_1 = hF(u)
       fft_2d=fft2(u);
       fft_1d=fft(u_one);
       %check the multiplications so that across rows is x, down columns is y
       totaltime=totaltime+dt;
       fft_u=fft2(u);
       xhat=real(ifft2(1i.*fft_u.*k'));
       yhat=real(ifft2(1i*k.*fft_u));
       q_1=dt*(-c.*xhat-d.*yhat);
       
       %second RK4 term q_2=hF(u+q_1/2)
       fft_u=fft2(u+q_1/2);
       xhat=real(ifft2(1i.*fft_u.*k'));
       yhat=real(ifft2(1i*k.*fft_u));
       q_2=dt*(-c.*xhat-d.*yhat);
       
       %third RK4 term q_3=hF(u+q_2/2)
       fft_u=fft2(u+q_2/2);
       xhat=real(ifft2(1i.*fft_u.*k'));
       yhat=real(ifft2(1i*k.*fft_u));
       q_3=dt*(-c.*xhat-d.*yhat);
       
       %fourth RK4 term q_4=hF(u+q_3)
       fft_u=fft2(u+q_3);
       xhat=real(ifft2(1i.*fft_u.*k'));
       yhat=real(ifft2(1i*k.*fft_u));
       q_4=dt*(-c.*xhat-d.*yhat);
              
       
       
       u_new=u+(q_1+2*q_2+2*q_3+q_4)/6;
       u=u_new;
       
       fft_u=fft(u_one);
       xhat=real(ifft(1i.*fft_u.*k_one));
       q_1=dt*(-c_one*xhat);
       
       %second RK4 term q_2=hF(u+q_1/2)
       fft_u=fft(u_one+q_1/2);
       xhat=real(ifft(1i.*fft_u.*k_one));
       q_2=dt*(-c_one*xhat);
       
       %third RK4 term q_3=hF(u+q_2/2)
       fft_u=fft(u_one+q_2/2);
       xhat=real(ifft(1i.*fft_u.*k_one));
       q_3=dt*(-c_one*xhat);
       
       %fourth RK4 term q_4=hF(u+q_3)
       fft_u=fft(u_one+q_3);
       xhat=real(ifft(1i.*fft_u.*k_one));
       q_4=dt*(-c_one*xhat);
              
       
       u_new_one=u_one+(q_1+2*q_2+2*q_3+q_4)/6;
       u_one=u_new_one;      
       if 0
       if mod(itr,10)==0
        surf(X,Y,u)
        axis([-pi+h pi -pi+h pi -2 2])
        %view(2)
        M(itr)=getframe;
       end
       end
end
error_inf(outer)=norm(u-init,inf);
toc
%clf
%surf(X,Y,u)
%axis([-pi+h pi -pi+h pi -0.1 2])
%view(2)
%hold on
%surf(X,Y,init)
end
if 0
loglog(Nvec,error_inf,'.')
hold on
semilogy(Nvec,Nvec.^(-4))
semilogy(Nvec,Nvec.^(-3))
end 