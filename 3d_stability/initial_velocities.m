%compute initial random velocities in Fourier space

function vel_i_hat = initial_velocities
	global N;
	%generate initial white noise solution
	for i=1:3
		vinit{i} = randn(N,N);
		vinit_hat{i}=modfft2(vinit{i});
	end

	for i=1:3
	    vel_i_hat{i}=vinit_hat{1}.*projection(1,i) + vinit_hat{2}.*projection(2,i) + vinit_hat{3}.*projection(3,i);	
	end

end
