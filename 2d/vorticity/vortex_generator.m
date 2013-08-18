clear all;

N=2^7;
h=2*pi/N;
x=-pi+h*(1:N); %linspace(-pi+h,pi,N);
y=x;   
[X,Y]=meshgrid(x,y);
omega=zeros(N);

for i=1:200
	x=-pi+rand*2*pi;
	y=-pi+rand*2*pi;
	i_x=15+rand*30;
	i_y=15+rand*30;

	if(rand>0.5)
		sign=1;
	else 
		sign=-1;
	end
	omega=omega+sign*exp(-i_x*(X+x).^2-i_y*(Y+y).^2);
end

surf(X,Y,omega,'EdgeColor','none')
axis([-pi+h pi -pi+h pi -10 10])
view(2);
dlmwrite('vortex',omega);
