function [angle_occupancy, theta] = DirectionalOccupancy(binsize,headdir,fs_video)
% Computes head direction occupancy
%
% Returns O, the directional occupancy of animal in Session root. If
% root.epoch is one epoch, O is a column vector. If root.epoch is more than
% one epoch, O is a matrix of column vectors of size M x N where M is the
% length of the number of bins as per binsize. N is size(root.epoch,1)
%
% Arguments:
% binsize -> positive number indicating binwidth in degrees
% Continuize -> 0 or 1 indicating whether to concatinate together multiple epochs
%
% Returns:
% O -> time in seconds in each bin
% theta -> bins in degrees
%
% O = root.DirectionalOccupancy;
% [O, theta] = root.DirectionalOccupancy(binsize, Continuize);

import CMBHOME.Utils.*

if ~exist('binsize', 'var'), binsize = 6; end

theta=-180+binsize/2:binsize:180-binsize/2;

theta = theta(:);

angle_occupancy = zeros(length(theta), size(headdir, 2), size(headdir, 3));

angle_occupancy = hist(headdir, theta) / fs_video ;
    
