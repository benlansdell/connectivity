function fig = plotmeanstd(fig, x, y, y_std, col, alph)
	%Note: to save transparency correctly use plot2svg instead of saveplot

	%Only accepts vector inputs currently
	if size(x,1) == 1
		x = x';
	end
	if size(y,1) == 1
		y = y';
	end
	if size(y_std,1) == 1
		y_std = y_std';
	end
	figure(fig);
	hold on
	%Make fill data to plot
	colwhite = col;%1-(1-col)*0.5;
	X = [x; flipud(x)];
	Y = [y+y_std; flipud(y-y_std)];
	vert = [X,Y];
	face = 1:size(vert,1);
	%Plot std area object as a whiter version of col
	%patch(X,Y,colwhite,'FaceColor', colwhite, 'FaceAlpha', alph);
	patch('Faces',face,'Vertices',vert,'FaceColor',colwhite, 'FaceAlpha', 0.1, 'EdgeAlpha', 0)
	%Plot mean curve
	%patchline(x, y, 'edgecolor', col, 'edgealpha', alph);
	plot(x, y, 'Color', col)
end