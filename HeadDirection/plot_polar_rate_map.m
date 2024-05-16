function [ratemap, ang_hd, mr, peak] = plot_polar_rate_map(Spks, Traj, msT, varargin)
    % root.plot_polar_rate_map(cel);
    % root.plot_polar_rate_map(cel, params, MaxRho);
    %
    % params is optional - [title, u2 score] (default [1 1])
    % MaxRho is optional - positive number, upper bound of Rho
    %
    % andrew 3 april 2010

    p = inputParser;
    p.addParamValue('fstart', 0, @isnumeric);
    p.addParamValue('animal', 'rat', @ischar);
    p.addParamValue('maxdist', 400, @isnumeric);
    p.addParamValue('khz', 25, @isnumeric);
    p.addParamValue('spon', 0, @isnumeric);
    p.addParamValue('speed', 2.5, @isnumeric);
    p.addParamValue('shuffle', 0, @isnumeric);
    p.addParamValue('shufflen', 0, @isnumeric);
    p.addParamValue('shuffletype', 'shift', @ischar);
    p.addParamValue('verbose', 0, @isnumeric);
    p.addParamValue('posture', [3 4 5 6], @isvector);
    p.addParamValue('dirtype', 'head', @ischar);

    p.parse(varargin{:});
    fstart = p.Results.fstart;
    animal = p.Results.animal;
    kHz = p.Results.khz;
    maxdist = p.Results.maxdist;
    spON = p.Results.spon;
    ThS = p.Results.speed;
    shuffle = p.Results.shuffle;
    shuffleN = p.Results.shufflen;
    shuffleType=p.Results.shuffletype;
    verbose = p.Results.verbose;
    posture = p.Results.posture;
    dirType = p.Results.dirtype;

    params = [1 1];
    MaxRho = 'yes';

    if size(Traj, 1) ~= size(msT, 1)

        if size(Traj, 1) < size(msT, 1)
            msT = msT(1:size(Traj, 1));
        else
            msFPS = median(diff(msT));
            seq = 1:length(Traj);
            seq = seq - 1;
            msT = seq * msFPS + msT(1);
            msT = msT';
        end

    end

    [ratemap, theta, wu2, ang_hd, mr, mvl, occumap, numSpks] = DirectionalTuningFnc(Spks, Traj, msT, 'fstart', fstart, 'animal', animal, 'kHz', kHz, 'spon', spON, 'speed', ThS, 'shuffle', shuffle,'shufflen',shuffleN,'shuffleType',shuffleType, 'posture', posture, 'maxdist', maxdist,'dirtype',dirType);

    theta = theta * unitsratio('rad', 'deg');
    theta = cat(1, theta(:), theta(1));
    theta = theta([length(theta) - floor(length(theta) / 4) + 1:length(theta) 1:length(theta) - floor(length(theta) / 4)]);

    if verbose == 1
        subplot(1, 3, 1);
        %plot_polar(ratemap,theta,params,MaxRho,wu2);
        plot_polar(ratemap, theta, params, MaxRho, wu2);
        subplot(1, 3, 2);
        plot_polar(numSpks, theta, params, MaxRho, wu2);
        subplot(1, 3, 3);
        plot_polar(occumap, theta, params, MaxRho, wu2);
    elseif verbose == 2
        plot_polar(ratemap, theta, params, MaxRho, wu2);
    end

    return;
    %%%%%%%%%%%%%
    function plot_polar(ratemap, theta, params, MaxRho, wu2)
        ratemap = cat(1, ratemap(:), ratemap(1));

        peak = max(ratemap);

        if ~exist('MaxRho', 'var')
            h1 = polar(theta(:), ratemap(:), 'k'); hold on;
        else

            if MaxRho == 'yes'
                %h1=polar(theta(:), ratemap(:), 'k'); hold on;
                h1 = myPolar(theta(:), ratemap(:), max(ratemap), 'k'); hold on;
            else
                h1 = myPolar(theta(:), ratemap(:), MaxRho, 'k'); hold on;
            end

        end

        set(h1, 'linewidth', 1.1)

        xs = ratemap(1:end - 1) .* cos(theta(1:end - 1)); % average
        ys = ratemap(1:end - 1) .* sin(theta(1:end - 1));

        coordlims = axis;

        ang_hd = atan2(mean(ys), mean(xs)); % mean direction

        mr = (cos(ang_hd) * sum(xs) + sin(ang_hd) * sum(ys)) / sum(ratemap(1:end - 1)); % mean resultant length

        mag_hd = sqrt(sum(ys) ^ 2 + sum(xs) ^ 2) / sqrt(sum(abs(ys)) ^ 2 + sum(abs(xs)) ^ 2) * coordlims(2); % for visualizations sake

        if ~exist('MaxRho', 'var')
            %    h1 = polar([ang_hd ang_hd], [0 mag_hd], 'r'); hold on
        else

            if MaxRho == 'yes'
                h1 = myPolar(theta(:), ratemap(:), max(ratemap), 'k'); hold on;
            else
                h1 = myPolar(theta(:), ratemap(:), MaxRho, 'k'); hold on;
            end

        end

        set(h1, 'linewidth', 1.1)

        xs = xlim;
        ys = ylim;

        if params(2)
            %        wu2 = HDWatsonsU2(); % calculate watsons u2
            if numel(wu2) == 1
                %            text(.55*xs(2), .8*ys(2), ['Watson''s U^2: ' num2str(wu2)], 'FontSize', 8);
            end

        end

        if params(1), title(['R:' num2str(mr, 3) ...
                                 'peak' num2str(peak, 3) 'Hz']); end

        return;
