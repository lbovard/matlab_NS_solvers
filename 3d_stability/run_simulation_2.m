clear all;

kz=[1.0:0.2:1.8];
Fh=0.2;
parfor i=1:length(kz);
	tic
	sigma(i)=simulation(kz(i),Fh);
	toc
end

dlmwrite('sigma_1.0_1.8',sigma)


