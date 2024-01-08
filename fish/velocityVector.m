function [firingRateMap, occupancyMap] = velocityVector(spikeCounts, directions, speeds, numBins, fs)
    debug = 0;
    % spikeCounts: Number of spikes at each time point
    % directions: Movement directions in degrees at each time point
    % speeds: Speed values at each time point
    % numXBins: Number of bins for x-direction speed
    % numYBins: Number of bins for y-direction speed
    % fs: Sampling frequency
    numXBins = numBins;
    numYBins = numBins;

    % Convert directions to radians
    directionsRad = deg2rad(directions);

    % Calculate X and Y components of speed
    xSpeeds = speeds .* cos(directionsRad);
    ySpeeds = speeds .* sin(directionsRad);

    % Initialize the bins for x and y direction speeds
    %xSpeedBins = linspace(min(xSpeeds), max(xSpeeds), numXBins + 1);
    %ySpeedBins = linspace(min(ySpeeds), max(ySpeeds), numYBins + 1);
  
    xSpeedBins = linspace(-100, 100, numXBins + 1);
    ySpeedBins = linspace(-100, 100, numYBins + 1);

    % Calculate occupancy map
    occupancyMap = histcounts2(xSpeeds, ySpeeds, xSpeedBins, ySpeedBins) / fs;

    % Initialize spike count map
    spikeCountMap = zeros(numXBins, numYBins);

    % Calculate spike counts for each x-y speed bin
    for i = 1:length(xSpeeds)
        xBin = find(xSpeedBins <= xSpeeds(i), 1, 'last');
        yBin = find(ySpeedBins <= ySpeeds(i), 1, 'last');

        if isempty(xBin)
            xBin = 1;
        end

        if isempty(yBin)
            yBin = 1;
        end

        if xBin > 0 && xBin <= numXBins && yBin > 0 && yBin <= numYBins
            spikeCountMap(xBin, yBin) = spikeCountMap(xBin, yBin) + spikeCounts(i);
        end

    end

    % Calculate firing rate map
    firingRateMap = spikeCountMap ./ occupancyMap;
    % Handle unvisited bins (set them to NaN)
    firingRateMap(occupancyMap == 0) = 0;

    firingRateMap = SmoothMat(firingRateMap, [9 9], 2.5);

    if debug
        % Display the firing rate map
        imagesc(firingRateMap);
        colormap('jet');
        colorbar;
        axis xy;
        title('Firing Rate in the Velocity-Vector Plane');
    end

end
