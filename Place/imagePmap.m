function imagePmap(rate_map, oc_map, max_rate)
    axis equal off;
    %ratemap

    if nargin < 3
        clims = [0 max(rate_map(:))];
        peak = max(rate_map(:));
    else
        clims = [0 max_rate];
        peak = max_rate;
    end

    no_occupancy = (oc_map == 0);
    [cbar, clims] = SmartColorbar(clims, 'jet(255)');
    rate_map(no_occupancy) = clims(1);

    imagesc(rate_map, clims);
    colormap(cbar);
    axis equal off

    gTitle = sprintf('Peak:%1.2f', peak);
    title(gTitle);

    return;
    %%%%%%%%%%%%%%%%%%%%%%%%
    function [cbar, clims] = SmartColorbar(clims, cmap_style)
        % finds the proper value for null data points to be plotted as white in an
        % IMAGESC plot, or similar.
        %
        % ARGUMENTS
        %   clims           (1x2 vector) either the range of the data, or the prefered range of the
        %                   colorbar
        %   cmap_style      string. ex. 'jet(255)'
        %
        % RETURNS
        %   cbar            Nx3 array of colorbar values. lowest is white
        %   clims           new clims. clims(1) should be included as the value of
        %                   the null data points to be properly mapped to white
        %
        % andrew 3 nov 2010
        % from CMBHOME

        eval(['cbar = ' cmap_style ';']); % initialize cbar

        cbar = cat(1, [1 1 1], cbar);

        nbins = size(cbar, 1) - 1;

        range = diff(clims) / (nbins - 1) + diff(clims);

        clims = [clims(2) - range, clims(2)];
