clear all;

kz=[[1:0.2:9.2], [9.4:0.1:10]];
Fh=0.2;
parfor i=1:length(kz);	
	sigma(i)=simulation(kz(i),Fh);
end

growth_data(1,:)=kz;
growth_data(2,:)=sigma;
dlmwrite('sigma_1_10.',sigma)
dlmwrite('sigma_1_10.full',growth_data);


