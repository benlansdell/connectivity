function y = softthreshold(X, thresh)
	y = sign(X).*max(0, abs(X)-thresh);
end