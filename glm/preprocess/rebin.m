function new_bs = rebin(bs, old_size, new_size)
	blocksize = new_size/old_size;
	assert(rem(blocksize,1) == 0, 'New bin size must be multiple of old binsize')
	n = floor(size(bs, 1)/blocksize);
	new_bs = zeros(n,1);
	for idx = 1:n
		new_bs(idx) = sum(bs((idx-1)*blocksize+1:(idx)*blocksize));
	end
end