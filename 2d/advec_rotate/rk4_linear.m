clear all;

Nvec=2.^(6);
for outer=1:length(Nvec)
N=Nvec(outer);
h=2*pi/N;
x=-pi+h*[1:N]; %linspace(-pi+h,pi,N);
%used for CLF condition
j=10;
dt=h/j;
tot=2*j*N;
totaltime=0;
%window=(tanh((X+pi-h)/0.1)-tanh((X-pi)/0.1)-1).*(tanh((Y+pi-h)/0.1)-tanh((Y-pi)/0.1)-1);
%advection co-efficients 
c=0.5;

%define the k matrix
k=[0:N/2-1 0 -N/2+1:-1];


tic
u=exp(-10*(x).^2);
%u=sin(X).*cos(Y);
init=u;

for itr=1:tot
       %first RK4 term q_1 = hF(u)
       
       %check the multiplications so that across rows is x, down columns is y
       totaltime=totaltime+dt;
       fft_u=fft(u);
       xhat=real(ifft(1i.*fft_u.*k));
       q_1=dt*(-c*xhat);
       
       %second RK4 term q_2=hF(u+q_1/2)
       fft_u=fft(u+q_1/2);
       xhat=real(ifft(1i.*fft_u.*k));
       q_2=dt*(-c*xhat);
       
       %third RK4 term q_3=hF(u+q_2/2)
       fft_u=fft(u+q_2/2);
       xhat=real(ifft(1i.*fft_u.*k));
       q_3=dt*(-c*xhat);
       
       %fourth RK4 term q_4=hF(u+q_3)
       fft_u=fft(u+q_3);
       xhat=real(ifft(1i.*fft_u.*k));
       q_4=dt*(-c*xhat);
              
       
       u_new=u+(q_1+2*q_2+2*q_3+q_4)/6;
       u=u_new;
       if 1
       if mod(itr,10)==0
        plot(x,u)
        axis([-pi+h pi 0 1])
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