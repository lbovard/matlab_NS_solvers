%precalculate the projection matrix

function precalc_projection
	global N k_x k_y k_z  k_inv_sq P
	p11=ones(N)-k_x.*k_x.*k_inv_sq;
    
	p12=-k_y.*k_x.*k_inv_sq;
	p13=-k_x.*k_z.*k_inv_sq;
	p22=ones(N)-k_y.*k_y.*k_inv_sq;
    
	p23=-k_y.*k_z.*k_inv_sq;
	p33=ones(N)-k_z.*k_z.*k_inv_sq;
    
	P=[[p11, p12, p13]; [p12, p22, p23]; [p13, p23, p33]];
end
