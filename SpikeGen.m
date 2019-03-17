

function [spikeMat] = SpikeGen(spikesPerS)

times = 500;	% a vector with each time step		
spikeMat = rand(1, times) < spikesPerS*0.1;
