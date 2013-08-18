%%
% This function accepts the fft of u,v (u_hat,v_hat) and the various wave numbers 
% computes f1 and f2 and p
%%
function [p_x,p_y,f1,f2,omega,p] = quantities(u,v,u_hat,v_hat,k_x,k_y,k_inv,k_sq)
	%obtain the various derivatives of u,v
	u_x=real(ifft2(1i.*k_x.*u_hat));
	u_y=real(ifft2(1i.*k_y.*u_hat));
	v_x=real(ifft2(1i.*k_x.*v_hat));
	v_y=real(ifft2(1i.*k_y.*v_hat));

	%compute the pressure
	f1=u.*u_x+v.*u_y;
	f2=u.*v_x+v.*v_y;
	f1hat=fft2(f1);
	f2hat=fft2(f2);
	p=real(ifft2(1i.*k_x.*f1hat.*k_inv+1i.*k_y.*f2hat.*k_inv));
	p_x = real(ifft2(1i.*k_x.*fft2(p)));
	p_y = real(ifft2(1i.*k_y.*fft2(p)));

	omega=v_x-u_y;
end
