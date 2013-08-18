function tru = modfft2(u)
	global N n_k cut;
	tru=fft2(u).*cut;
end
