N=128;
x=h*(1:N);
h=2*pi/N;
t=0;
dt=h/4;
c=.2+sin(x-1).^2;
v=exp(-100*(x-1).^2); vold=exp(-100*(x-.2*dt-1).^2);

tmax=20;
tplot=0.15;
clf, drawnow
plotgap=round(tplot/dt); dt=tplot/plotgap;
nplots=round(tmax/tplot);
data=[v;zeros(nplots,N)]; tdata=t;
for i=1:nplots
    for n=1:plotgap
        t=t+dt;
        v_hat=fft(v);
        w_hat=1i*[0:N/2-1 0 -N/2+1:-1].*v_hat;
        w=real(ifft(w_hat));
        vnew=vold-2*dt*c.*w; vold=v; v=vnew;
    end
    data(i+1,:)=v; tdata=[tdata;t];
end
waterfall(x,tdata,data), view(10,70),colormap([0 0 0])
axis([0 2*pi 0 tmax 0 5]), ylabel t, zlabel u, grid off