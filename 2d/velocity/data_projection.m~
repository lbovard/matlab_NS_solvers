function [F_hat,G_hat,omega]=data_projection(u,v,u_hat,v_hat,k_x,k_y,k_inv,k_sq,P_11,P_12,P_22)

	%compute derivatives
	u_x=real(ifft2(1i*k_x.*u_hat));
	u_y=real(ifft2(1i*k_y.*u_hat));
	v_x=real(ifft2(1i*k_x.*v_hat));
	v_y=real(ifft2(1i*k_y.*v_hat));
	f1=u.*u_x+v.*u_y;
	f2=u.*v_x+v.*v_y;
	f1_hat=fft2(f1);
	f2_hat=fft2(f2);
	
	%symm of projection tensor
	P_21=P_12;
	%update 
	F_hat=-P_11.*f1_hat-P_12.*f2_hat;
	G_hat=-P_21.*f1_hat-P_22.*f2_hat;

	omega=v_x-u_y;
end
