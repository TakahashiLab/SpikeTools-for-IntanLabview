%tagPlot
%significant interneuron or pyramidal cells must met two criteria described below.
%1. The firing rate during control period (presilent period) was 1 Hz or more.
%2. The firing rate during optogenetics stimulation period  was above the control.
function [c, raster] = tagPlot(unit, tagPoint)

    kHz = 25;
    dispRangePre = kHz * 200; %ms before stimOn
    dispRangePost = kHz * 200; %ms after stimOff
    stimDuration = kHz * 50; %ms
    %stimOff=dispRangePre+kHz*20;
    binWidth = 10; %ms
    step = binWidth / (1 / kHz);
    alpha = 0.05;
    frTh = 1.0; % 1.0Hz or more

    triggerPoint = tagPoint;

    clf;
    raster = [];
    Ind = [];
    cnt = 1;
    rasterOrg = [];

    tmps = unit;
    tmps = double(tmps);
    tmpsMs = floor(tmps / kHz);

    for i = 1:length(triggerPoint)
        raster = [raster tmps - triggerPoint(i)];
        Ind = [Ind ones(size(tmps)) * cnt];

        rasterOrg = [rasterOrg tmpsMs - floor((triggerPoint(i) - dispRangePre)/kHz)];

        cnt = cnt + 1;
    end

    IndSp = Ind;
    IndSp(rasterOrg <= 0) = [];
    rasterOrg(rasterOrg <= 0) = [];
  
    spt = sparse(IndSp,rasterOrg,  1);

    stimAfter = dispRangePre + stimDuration;

    subplot(2, 1, 1);
    plot(raster, Ind, '.');
    axis([-dispRangePre stimDuration + dispRangePost 0 cnt]);
    tics = [-dispRangePre 0 stimDuration stimDuration + dispRangePost];
    hold on;
    plot([0 0], [0 cnt], 'r');
    plot([stimDuration stimDuration], [0 cnt], 'g');
    set(gca, 'xtick', tics, 'xticklabel', tics / kHz);

    subplot(2, 1, 2);
    %raster(raster <- dispRangePre) = [];
    %raster(raster > (stimDuration + dispRangePost)) = [];

    for i = 1:max(Ind)

        q = histcounts(raster(Ind == i), -dispRangePre:step:(stimDuration + dispRangePost));
        q(q > 1) = 1;

        if i == 1
            c = q;
        else
            c = c + q;
        end

    end

    cOrg = c;
    c = (c ./ (cnt - 1)) ./ (binWidth / 1000);
    bar(c);

    bf = mean(c(1:dispRangePre / step));
    sf = mean(c(dispRangePre / step:stimAfter / step));

    fprintf('base fr=%f, stim fr=%f\n', bf, sf);

    %poisson test
    c = double(c);

    [pInt, pPyr, pPyrReb] = poissonTest(cOrg((dispRangePre / step) + 1:stimAfter / step), cOrg(1:(dispRangePre / step) - 1));
 
   
    dt = 1 ; %ms
    win = [-200 50];
    [p_value, Idiff] = isikldist(spt, dt, win, binWidth);

    pInt=p_value;
    if bf > frTh

        if pInt < alpha
            fprintf('Int:: isikldist test=%f\n', pInt);
        elseif pPyr < alpha
            fprintf('Pyr:: poisson test=%f\n', pPyr);
        elseif pPyrReb < alpha
            fprintf('PyrReb:: poisson test=%f\n', pPyrReb);
        end

    end

    maxC = max(c);

    if maxC == 0
        maxC = 1;
    end

    axis([0 (dispRangePre + stimDuration + dispRangePost) / step 0 maxC * 1.2]);
    hold on;

    plot([dispRangePre / step dispRangePre / step], [0 max(c) * 1.2], 'r');
    plot([stimAfter / step stimAfter / step], [0 max(c) * 1.2], 'g');

    tics = [0 dispRangePre / step stimAfter / step];
    ticlabel = [-dispRangePre / step * binWidth 0 stimDuration / step * binWidth];
    set(gca, 'xtick', tics, 'xticklabel', ticlabel);

    return;
    %%%%%%%%%
    function [pInt, pPyr, pPyrReb] = poissonTest(X, preX)
        mu = mean(preX);

        pInt = poisscdf(X(1), mu, 'upper');
        pPyr = poisscdf(X(3:end), mu) .* (length(X) - 2);
        pPyr = min(pPyr);
        pPyrReb = poisscdf(X(2:end), mu, 'upper') .* (length(X) - 1);
        pPyrReb = min(pPyrReb); %rebound
        return;
        %%%%%%%%%%%%
        kHz = 25;
        window = 500; %msec
        win = window * kHz;

        loop = size(tagPoint, 2);
        unit = unit(unit > tagPoint(1) - win & unit < tagPoint(end) + win);
        Spks = unit;
        trials = unit;

        for i = 1:loop
            spk = unit - tagPoint(i);
            ind = find(spk >- win & spk < win);
            Spks(ind) = unit(ind) - tagPoint(i);
            trials(ind) = i;
        end

        plot(Spks, trials, '.');
        return;
