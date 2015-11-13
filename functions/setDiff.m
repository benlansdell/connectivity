function difference = setDiff(A, B)
	anotinb = A(~ismember(A,B));
	bnotina = B(~ismember(B,A));
	difference = sort([anotinb, bnotina]);
end