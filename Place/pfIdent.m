function [PeakPos] = pfIdent(rm, oc, Traj, spatial_scale, binside, verbose)
    ThF = 0;
    binside2 = binside ^ 2;
    thRate = 0.5;

    sp_thr = max(max(rm)) / 3; %not 3

    b = bwboundaries(oc);
    b = b{1};

    rm_bord = zeros(size(rm, 1), size(rm, 2));

    for i = 1:length(b)
        rm_bord(b(i, 1), b(i, 2)) = 1;
    end

    field = (rm >= sp_thr);

    a = regionprops(field, {'Area' 'PixelList'});
    b = [];

    c = 1;

    for i = 1:length(a)

        if ((a(i).Area * binside2) > ThF)
            b{c} = a(i).PixelList;
            c = c + 1;
        end

    end

    if verbose
        subplot(2, 2, 1);
        imagePmap(rm, oc);
        axis on;
    end

    offsetX = min(Traj(:, 3));
    offsetY = min(Traj(:, 4));
    offSets = [offsetX offsetY];
    PeakPos = [];
    orgPeakPos = [];

    for i = 1:c - 1
        xy = b{i};
        ind = sub2ind(size(rm), xy(:, 2), xy(:, 1));
        [maxRate, ind] = max(rm(ind));

        if maxRate > thRate
            orgPeakPos = [orgPeakPos; xy(ind, :)];

            pxy = xy;
            xy = xy .* (1 / spatial_scale) .* binside;
            xy(:, 1) = xy(:, 1) + offsetX;
            xy(:, 2) = xy(:, 2) + offsetY;

            PeakPos = [PeakPos; pxy(ind, :)];

            if verbose
                subplot(2, 2, 1);
                hold on;
                plot(pxy(ind, 1), pxy(ind, 2), 'k*');
                axis equal;

            end

        end

    end

    return;
