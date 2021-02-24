function [rate_map, xdim, occupancy, no_occupancy] = rateLinearMap(spk_x,x,spatial_scale,fs_video,fieldShuffleON)

binside=2.5;
GRID=[9 9];

std_smooth_kernel=5;
%[occupancy, xdim, ydim] = Occupancy(x,y,
%spatial_scale,binside,fs_video);
dim=binside;;

xdim = min(x):spatial_scale^-1*dim:max(x); %edges of x and y
                                           %dimensions
occupancy = histcounts(x, xdim)/fs_video;
no_occupancy = occupancy==0; 
%spikes = hist3([spk_x, spk_y], 'Edges', {xdim, ydim});
spikes = histcounts(spk_x, xdim);
rate_map=spikes./occupancy;
rate_map(no_occupancy) = 0; % set no occupancy to zero
rate_map=SmoothMat(rate_map,GRID,binside);
%rate_map = SmoothMat(spikes, 10, 2)./SmoothMat(occupancy, 10, 2); 

rate_map = rate_map'; % returns these three
occupancy = occupancy';
no_occupancy = no_occupancy';


return;

%%%%%%%%%%%%%%%%%%%%%%%%%5555555555
function mat = SmoothMat(mat, kernel_size, std)
if nargin<3
    std=1;
end

if std == 0, return; end


Rgrid=-kernel_size/2: kernel_size/2;
kernel=normpdf(Rgrid,0,std);

kernel = kernel./sum(sum(kernel));

mat = conv(mat, kernel, 'same');
return;
