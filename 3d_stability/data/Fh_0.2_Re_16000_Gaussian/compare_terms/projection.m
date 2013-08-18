%returns the projection tensor
%replace with cell array:
function p = projection(i,j)
	global P N
	%selects correct values in giant P matrix
	a=(i-1)*N+1:i*N;
	b=(j-1)*N+1:j*N;
	p=P(a,b);
end
