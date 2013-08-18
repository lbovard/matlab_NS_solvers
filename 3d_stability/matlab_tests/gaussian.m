% Function that initialises lamb dipole

%fix up velocities to remove redundant calculations
function gaussian 
        global X Y N u_0 v_0 om_i alpha k_x k_y;
	k_sq=k_y.^2+k_x.^2;
	kzeros=find(k_sq==0);
	k_inv=1./k_sq;
	k_inv(kzeros)=0;
        omega=exp(-(X-1/(2*alpha)).^2-Y.^2)+exp(-(X+1/(2*alpha)).^2-Y.^2);
	om_i=omega;
	fft_omega=modfft2(omega);
	u_0=ifft2(1i*k_y.*fft_omega.*k_inv);
	v_0=ifft2(-1i*k_x.*fft_omega.*k_inv);
	surf(X,Y,real(omega));
end
