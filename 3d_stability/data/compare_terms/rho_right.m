%compute right side of rho function

function rr= rho_right(w_hat,P_hat,t)
	global Fh u_0 v_0 k_x k_y k_sq Re Sc;
	
	rho=ifft2(exp(-k_sq*t/Re/Sc).*P_hat);

	rr=exp(k_sq*t/Re/Sc).*(-1i*k_x.*modfft2(u_0.*rho)-1i*k_y.*modfft2(v_0.*rho)+w_hat/Fh^2);
end
