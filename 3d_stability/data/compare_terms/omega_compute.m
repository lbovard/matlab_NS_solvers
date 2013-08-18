%compute the vorticity and convert to real space

function omega = omega_compute(vel_hat)
	global k_x k_y k_z 
	omega{1} =ifft2(1i*k_y.*vel_hat{3}-1i*k_z.*vel_hat{2});
	omega{2} =ifft2(1i*k_z.*vel_hat{1}-1i*k_x.*vel_hat{3});
	omega{3} =ifft2(1i*k_x.*vel_hat{2}-1i*k_y.*vel_hat{1});
end
