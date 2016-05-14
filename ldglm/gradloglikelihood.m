%These assume an exponential non-linearity
function G = gradloglikelihood(S, Y)
	%Luckily this can be done element-wise since each s component
	%is conditionally independent of the others given y.
	G = S - exp(Y);
end