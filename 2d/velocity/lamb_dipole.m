% Function that initialises lamb dipole

%fix up velocities to remove redundant calculations

N=2^7;

h=2*pi/N;
x=-pi+h*(1:N); %linspace(-pi+h,pi,N);
y=x;
[X,Y]=meshgrid(x,y);
%define the k matrix
k_x=repmat([0:N/2-1 0 -N/2+1:-1],N,1);
k_y=k_x';
k_sq=k_y.^2+k_x.^2;
k_inv=1./k_sq;
%hackish method to find the div by 0, replace
k_zeros=find(k_inv==Inf);
k_inv(k_zeros)=0;

%zero of J_{1}(b)=0, more accuracy?
r=sqrt(X.^2+Y.^2);
spec_value=find(r==0);
r_sq=r.^2;
theta=atan2(Y,X);
b=3.8317;
%need J_{1}'(b), use J_{1}'(b) = J_{0}(b) (see 59 of mathworld entry of
%Bessel functions of the first kind)
A=besselj(0,b);

omega=zeros(N,N);       
omega=2*b/A*besselj(1,b*r).*sin(theta); 
k =find(r>1);
omega(k)=0;
om_i=omega;
fft_omega=fft2(omega);

u_0=zeros(N,N);
v_0=zeros(N,N);


u_0=2./(A*b*r.^3).*((X.^2-Y.^2).*besselj(1,b*r)+b*Y.^2.*r.*besselj(0,b*r));
v_0=-2.*X.*Y./(r.^3*b*A).*(-2*besselj(1,b*r)+b.*r.*besselj(0,b*r));
%matlab cannot compute the limit
u_0(spec_value)=1/(A);
v_0(spec_value)=0;


u_0(k)=(1-1./r_sq(k))+2*Y(k).^2./r_sq(k).^2;
v_0(k)=-2*X(k).*Y(k)./r_sq(k).^2;
v=-1i*k_x.*fft_omega.*k_inv;
v(1,1)=mean(v_0(:));
v=real(ifft2(v));
u=1i*k_y.*fft_omega.*k_inv;
%u(1,1)=mean(u_0(:))*N^2;
u=real(ifft2(u))+0.9204*ones(N,N);

dlmwrite('dipolev',omega)
