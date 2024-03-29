function [rate_map, xdim, ydim, occupancy, no_occupancy] = ratemap(spk_x, spk_y, x, y, spatial_scale, fs_video, binside, GRID, xdim, ydim)

    dim = 1;
    omit = 1;

    if nargin == 6
        binside = 2.5;
        GRID = [9 9];
        [occupancy, xdim, ydim] = Occupancy(x, y, spatial_scale, dim, fs_video);
    elseif nargin == 7
        GRID = [9 9];
        [occupancy, xdim, ydim] = Occupancy(x, y, spatial_scale, dim, ...
            fs_video);
    elseif nargin == 8

        if isnumeric(GRID)
            GRID = [9 9];
            omit = 0;
        end

        [occupancy, xdim, ydim] = Occupancy(x, y, spatial_scale, dim, ...
            fs_video);
    else
        [occupancy] = Occupancy(x, y, spatial_scale, dim, fs_video, xdim, ydim);
    end

    if omit
        occupancy = OmitIslands(occupancy);
    end

    no_occupancy = occupancy == 0;
    spikes = hist3([spk_x, spk_y], 'Edges', {xdim, ydim});

    rate_map = spikes ./ occupancy;
    rate_map(no_occupancy) = 0; % set no occupancy to zero
    rate_map = SmoothMat(rate_map, GRID, binside);
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
        if nargin < 3
            std = 1;
        end

        if std == 0, return; end

        [Xgrid, Ygrid] = meshgrid(-kernel_size(1) / 2:kernel_size(1) / 2, -kernel_size(2) / 2:kernel_size(2) / 2);

        Rgrid = sqrt((Xgrid .^ 2 + Ygrid .^ 2));

        kernel = pdf('Normal', Rgrid, 0, std);

        kernel = kernel ./ sum(sum(kernel));

        mat = conv2(mat, kernel, 'same');
        return;
