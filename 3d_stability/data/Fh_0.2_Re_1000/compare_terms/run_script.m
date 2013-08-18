% 1 = u_diffusion, 2=v_diffusion, 3=w_diffusion, 4=rho_diffusion
% 5 = u_current_advec, 6 = v_current_advec, 7 = w_currect_advec
% 8 = u_basic_advec, 9 = v_basic_advec, 10= w_basic_advec
% 11= u_rho_proj, 12 = v_rho_proj, 13=w_rho_proj
% 14 = rho_u, 15 = rho_v, 16=rho_w

clear all;
clf;
close all;
kz=[1.4:2:10, 11:1:15, 20:5:50];
Fh=0.2;
get_values=0;
plots=0;
if(get_values==1)
for i=1:length(kz)
	term_norms(i,1:16)=compare_terms(kz(i),Fh,plots);
	if(mod(i,10)==0) 
		display(i) 
	end
end
else
	term_norms=dlmread('term_norms');
end
dlmwrite('term_norms',term_norms);

u_data=term_norms(:,[1,5,8,11]);
v_data=term_norms(:,[2,6,9,12]);
w_data=term_norms(:,[3,7,10,13]);
rho_data=term_norms(:,[4,14,15,16]);

data{1}=u_data;
data{2}=v_data;
data{3}=w_data;
data{4}=rho_data;
for j=1:4
for i=1:length(kz)
	[a,ind]=max(data{j}(i,:));
	max_term(i)=ind;
	data{j}(i,:)=data{j}(i,:)/data{j}(i,max_term(i));
end
end
figure 
eq1=subplot(2,2,1);
plot(kz,data{1}) 
title('u equation')
xlabel('k_z')
legend('diff','cadvec', 'badvec', 'rho','Location','East');

eq2=subplot(2,2,2);
plot(kz,data{2})
title('v equation')
xlabel('k_z')
legend('diff','cadvec', 'badvec', 'rho','Location','East');

eq3=subplot(2,2,3);
plot(kz,data{3})
title('w equation')
xlabel('k_z')
legend('diff','cadvec', 'badvec', 'rho','Location','East');
eq4=subplot(2,2,4);
plot(kz,data{4})
title('rho equation')
xlabel('k_z')
legend('diff','rho u', 'rho v', 'rho w','Location','East');


linkaxes([eq1 eq2 eq3 eq4], 'xy')
