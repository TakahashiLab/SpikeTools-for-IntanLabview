%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [occupancy, xdim, ydim] = Occupancy(x,y,spatial_scale,dim,fs_video)
% Returns occupancy matrix of time spent in 2D spatial bins 
%
% Returns occupancy matrix of time spent in bins of dimension 3x3cm^2 be
% default. occupancy = root.Occupancy(xdim, ydim) where will return a matrix of
% dimensions xdim by ydim(:,2)
%
% Remember that root.spatial_scale indicates cm/pixel in your experiment.
%
% Occupancy is only returned for data within root.epoch. If multiple epochs
% are selected, they are all concatinated and one matrix is returned.
% [occupancy, xdim, ydim] = root.Occupancy returns occupancy matrix and x
% and y dimensions
%
% andrew bogaard 17 may 2010
% from https://github.com/hasselmonians/CMBHOME
xdim = min(x):spatial_scale^-1*dim:max(x); %edges of x and y dimensions
ydim = min(y):spatial_scale^-1*dim:max(y);

occupancy = hist3([x, y], 'Edges', {xdim, ydim}) /fs_video;
return;