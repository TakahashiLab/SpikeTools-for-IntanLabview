function [rate_map, xdim, ydim, occupancy, no_occupancy] = ratemap(spk_x,spk_y,x,y,spatial_scale,fs_video,fieldShuffleON)

binside=2.5;

GRID=[9 9];
[occupancy, xdim, ydim] = Occupancy(x,y, spatial_scale,binside,fs_video);
occupancy = OmitIslands(occupancy);

no_occupancy = occupancy==0; 
spikes = hist3([spk_x, spk_y], 'Edges', {xdim, ydim});

%rate_map = SmoothMat(spikes, [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside)./SmoothMat(occupancy, [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside); 

rate_map=spikes./occupancy;
rate_map(no_occupancy) = 0; % set no occupancy to zero
rate_map=SmoothMat(rate_map, GRID,binside);
rate_map = rate_map'; % returns these three
occupancy = occupancy';
no_occupancy = no_occupancy';
    
return;

%%%%%%%%%%%%%%%%%%%%%%%%%5555555555
function mat = SmoothMat(mat, kernel_size, std)
%
% Smooths matrix by convolving with 2d gaussian of size
% kernel_size=[bins_x bins_y] and standard deviation 'std'
%
% if std==0, just returns mat
%
% 10 december 2009 andrew
%
% from https://github.com/hasselmonians/CMBHOME
if nargin<3
    std=1;
end

if std == 0, return; end

[Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);

Rgrid=sqrt((Xgrid.^2+Ygrid.^2));

kernel = pdf('Normal', Rgrid, 0, std);

kernel = kernel./sum(sum(kernel));

mat = conv2(mat, kernel, 'same');
return;
