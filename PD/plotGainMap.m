function outMph = plotGainMap(phaseHist, varargin)

    p = inputParser;
    p.addParamValue('display', 'both', @ischar);
    p.addParamValue('band', 'all', @ischar);
    p.addParamValue('hist', 'off', @ischar);
    p.addParamValue('topfreq', 120, @isnumeric);

    p.parse(varargin{:});

    disp = p.Results.display;
    band = p.Results.band;
    histon = p.Results.hist;
    topFreq = p.Results.topfreq;

    sd = 2; %5
    edges = [-3 * sd:1:3 * sd];
    kernel = normpdf(edges, 0, sd) * sd;
    kernel = kernel ./ sum(kernel(:));

    %imagesc(conv2(kernel,kernel,mean(phaseHist(:,[1:40 80:120],:),3)','same'));
    %mph=conv2(kernel,kernel,mean(phaseHist(:,[1:40 81:120],:),3)','same');

    switch band
        case 'theta',
            frange = 4:11;
        case 'beta',
            frange = 12:30;
        case 'gamma',
            frange = 31:80;
            otherwise,
            frange = [];
    end

    %topFreq=80;
    %topFreq = 120;
    topFreq = topFreq / 2;

    switch lower(disp)
        case 'both',

            %mph=conv2(kernel,kernel,mean(phaseHist(:,[1:40],:),3)','same');
            %mph=SmoothMat(mean(phaseHist(:,[1:40],:),3)', [9 9],2);
            %mph=SmoothMat(mean(phaseHist(:,[1:topFreq],:),3)', [9 9],2);
            mph = mean(phaseHist(:, [1:topFreq], :), 3)';
            mph = [mph(:, 11:20) mph mph(:, 1:10)];
            mph = SmoothMat(mph, [9 9], 2);
            subplot(1, 2, 1);
            %imagesc([mph(:,11:20) mph mph(:,1:10)]);
            imagesc(mph);
            title('ascending');
            set(gca, 'ydir', 'normal');
            set(gca, 'ytick', [1:5:topFreq], 'yticklabel', [0:5:topFreq] * 2);
            set(gca, 'xtick', [1 11 21 31 41], 'xticklabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});
            hold on;
            t = 0:40;
            plot(t + 1, cos(t / pi) * topFreq / 2 + topFreq / 2, 'w');
            plot([0 40], [4/2 4/2], 'k');
            plot([0 40], [12/2 12/2], 'k');
            plot([0 40], [30/2 30/2], 'k');
            plot([0 40], [80/2 80/2], 'k');
            plot([21 21], [0 topFreq], 'k');

            colorbar;
            xlabel('Phase [radian]');
            ylabel('Frequency [Hz]');

            subplot(1, 2, 2);
            %%%mph=conv2(kernel,kernel,mean(phaseHist(:,[81:120],:),3)','same');
            %mph=conv2(kernel,kernel,mean(phaseHist(:,[120:-1:81],:),3)','same');
            %mph=SmoothMat(mean(phaseHist(:,[120:-1:81],:),3)', [9 9],2);
            %mph=SmoothMat(mean(phaseHist(:,[120:-1:120-topFreq+1],:),3)', [9 9],2);
            mph = mean(phaseHist(:, [120:-1:120 - topFreq + 1], :), 3)';
            mph = SmoothMat([mph(:, 11:20) mph mph(:, 1:10)], [9 9], 2);
            %    imagesc([mph(:,11:20) mph mph(:,1:10)]);
            imagesc(mph);
            title('descending');
            set(gca, 'ydir', 'normal');
            %set(gca,'ytick',[1:5:40],'yticklabel',[(40:-5:1)-4]*2);
            set(gca, 'ytick', [1:5:topFreq], 'yticklabel', [0:5:topFreq] * 2);
            set(gca, 'xtick', [1 11 21 31 41], 'xticklabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});
            hold on;
            t = 0:40;
            plot(t + 1, cos(t / pi) * topFreq / 2 + topFreq / 2, 'w');
            xlabel('Phase [radian]');
            ylabel('Frequency [Hz]');

            plot([0 40], [4/2 4/2], 'k');
            plot([0 40], [12/2 12/2], 'k');
            plot([0 40], [30/2 30/2], 'k');
            plot([0 40], [80/2 80/2], 'k');
            plot([21 21], [0 topFreq], 'k');

            colormap(jet);
            colorbar;
            caxis([0.5 max(mph(:))]);

        case 'ascend',
            xr = [1:topFreq];
            xrange = [0:5:topFreq];
            outMph = plotfnc(phaseHist, xr, xrange, kernel, frange, band, histon);

        case 'descend',
            xr = [120:-1:120 - topFreq + 1]; %81
            xrange = [0:5:topFreq];

            outMph = plotfnc(phaseHist, xr, xrange, kernel, frange, band, histon);
            %%title('descending');
    end

    return;
    %%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function outMph = plotfnc(phaseHist, xr, xrange, kernel, frange, band, histon)
    binsize = 2;
    %topFreq = 80;
    %topFreq = topFreq / 2;
    %xr = [1:topFreq];
    %xrange = [0:5:topFreq];
    outMph = [];

    if strcmp(band, 'all')
        mph = mean(phaseHist(:, xr, :), 3)';
        mph = [mph(:, 11:20) mph mph(:, 1:10)];
    elseif strcmp(band, 'phase')
        mph = squeeze(mean(phaseHist(:, xr, :), 2));
        outMph = mph;
        seMph = (std(mph, 0, 2) ./ sqrt(size(mph, 2))')'; %sem
        mph = mean(mph, 2)';

        mph = mph(1:end - 1);
        seMph = seMph(1:end - 1);
        outMph = outMph(1:end - 1, :);

        hf = length(mph) / 2;
        mph = [mph(:, hf + 1:end) mph mph(:, 1:hf)];
        seMph = [seMph(:, hf + 1:end) seMph seMph(:, 1:hf)];
        [~, ~, bins] = histcounts(1:length(mph), length(mph) / binsize);
        mph = accumarray(bins', mph, [], @mean)';
        seMph = accumarray(bins', seMph, [], @mean)';

        [~, ~, bins] = histcounts(1:size(outMph, 1), size(outMph, 1) / binsize);

        binOutMph = zeros(size(outMph));
        binOutMph(size(binOutMph, 1) / binsize + 1:end, :) = [];

        for q = 1:size(outMph, 2)
            binOutMph(:, q) = accumarray(bins', outMph(:, q), [], @mean)';
        end

        outMph = binOutMph;

        plot(mph);

        hold on; plot(mph + seMph, 'r'); plot(mph - seMph, 'r');

        set(gca, 'xtick', [0:10 / binsize:40 / binsize], 'xticklabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});
        xlabel('phase');
        ylabel('Gain');

    elseif strcmp(band, 'freq')

        mph = squeeze(mean(phaseHist(:, xr, :), 1));
        outMph = mph;
        mph = mean(mph, 2)';
        [~, ~, bins] = histcounts(1:length(mph), length(mph) / binsize);
        mph = accumarray(bins', mph, [], @mean)';
        %mph=SmoothMat(mph,[3 3],1);
        plot(mph);

        set(gca, 'xtick', xrange, 'xticklabel', xrange * 2 * binsize);
        xlabel('Frequency [Hz]');
        ylabel('Gain');
    else
        frange = unique(floor(frange / 2));

        mph = squeeze(mean(phaseHist(:, xr(frange), :), 2))';
        [maxMph, maxMphId] = max(mph, [], 2);
        [~, id] = sort(maxMphId, 'descend');
        %mph = mph(:, id);

        [~, id] = sort(maxMph, 'descend');
        mph = mph(id,:);
        mph = [mph(:, 11:20) mph mph(:, 1:10)];
    end

    if (strcmp(band, 'freq') | strcmp(band, 'phase'))
     
    elseif strcmp(histon, 'off')

        mph = SmoothMat(mph, [9 9], 2);
        %mph=conv2(kernel,kernel,mph,'same');
        imagesc(mph);

        %%title('ascending');
        set(gca, 'ydir', 'normal');

        if strcmp(band, 'all')
            set(gca, 'ytick', xrange, 'yticklabel', xrange * 2);
            ylabel('Frequency [Hz]');
        else

            if size(mph, 1) ~= 1
                set(gca, 'ytick', [1 size(mph, 1)], 'yticklabel', [1 size(mph, 1)]);
            end

            ylabel('Neuron number');
        end

        set(gca, 'xtick', [1 11 21 31 41], 'xticklabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});
        hold on;
        t = 0:40;
        plot(t + 1, cos(t / pi) * 10 + 10.5, 'w');
        xlabel('Phase [radian]');
        colormap(jet);
        colorbar;

        if size(mph, 1) ~= 1
            caxis([0.5 max(mph(:))]);
        end

    else
        mph2 = [];
        bin = 1;
        sd = 1;
        edges = [-3 * sd:3 * sd];
        kernel = normpdf(edges, 0, sd) * sd;
        %mph=conv(mean(mph),kernel,'same');
        mph = mean(mph);
        mph = SmoothMat(mph, [3 3], 1);
        bar(mph);
        set(gca, 'xtick', [1 11 21 31 41], 'xticklabel', {'-2\pi', '-\pi', '0', '\pi', '2\pi'});
    end

    return;

end
