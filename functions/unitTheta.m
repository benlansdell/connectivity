function [direction, torquesize] = unitTheta(torqueRU, torqueFE)	
	%Returns angle of 2D vector (or list of 2D vectors)
	torquesize = sqrt(torqueRU.^2+torqueFE.^2);
	%Compute angle
	torqueH = atan(torqueRU./torqueFE); %theta
	torquesize = sqrt(torqueRU.^2 + torqueFE.^2); %speed
	piplus = torqueFE < 0 & torqueRU >= 0;
	piminus = torqueFE < 0 & torqueRU < 0;
	%Between plus/minus pi
	torqueH(piplus) = torqueH(piplus) + pi;
	torqueH(piminus) = torqueH(piminus) - pi;
	%Between 0 and 2pi
	torqueH(torqueH<0) = torqueH(torqueH<0)+2*pi;
	direction = torqueH;
end