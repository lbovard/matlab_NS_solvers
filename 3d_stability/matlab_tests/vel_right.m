%compute right hand side of velocity equation

function vr=vel_right(U_hat,P_hat,t)
	global om_i u_0 v_0 k_sq Re Sc;

	for i=1:3
		vel{i}=ifft2(U_hat{i}.*exp(-k_sq*t/Re));
		vel_hat{i}=modfft2(vel{i});
	end

	rho_hat=P_hat.*exp(-k_sq*t/Re/Sc);

	omega = omega_compute(vel_hat);
	A=vel{2}.*om_i+v_0.*omega{3};
	B=-vel{1}.*om_i-u_0.*omega{3};
	C=u_0.*omega{2}-v_0.*omega{1};
	MM{1} = modfft2(A);
	MM{2} = modfft2(B);
	MM{3} = modfft2(C) - rho_hat;

	for i=1:3
		vr{i} = exp(k_sq*t/Re).*(projection(i,1).*MM{1} + projection(i,2).*MM{2} + projection(i,3).*MM{3});
	end			
end
