function [F,G,omega]=data(u,v,u_hat,v_hat,k_x,k_y,k_inv,k_sq)

	%compute derivatives
	u_x=real(ifft2(1i*k_x.*u_hat));
	u_y=real(ifft2(1i*k_y.*u_hat));
	v_x=real(ifft2(1i*k_x.*v_hat));
	v_y=real(ifft2(1i*k_y.*v_hat));
	f1=u.*u_x+v.*u_y;
	f2=u.*v_x+v.*v_y;
	f1_hat=fft2(f1);
	f2_hat=fft2(f2);
	p_hat=1i*k_x.*f1_hat.*k_inv+1i*k_y.*f2_hat.*k_inv;
	p_x=real(ifft2(1i*k_x.*p_hat));
	p_y=real(ifft2(1i*k_y.*p_hat));
	F=fft2(-f1-p_x);
	G=fft2(-f2-p_y);
	
	omega=v_x-u_y;
end
