%
function [normalPyrS, normalIntS, pdPyrS, pdIntS] = Gainall(basename, varargin)
    p = inputParser;
    p.addParamValue('proc', 'individual', @ischar);
    p.addParamValue('stimulation', 'both', @ischar);
    p.addParamValue('display', 'ascend', @ischar);
    p.addParamValue('plot', 'all', @ischar);
    p.addParamValue('band', 'non', @ischar);
    p.addParamValue('hist', 'off', @ischar);
    p.addParamValue('topfreq', 120, @isnumeric);
    p.addParamValue('alpha', 0.05, @isnumeric);

    p.parse(varargin{:});

    proc = p.Results.proc;
    stim = p.Results.stimulation;
    disp = p.Results.display;
    plt = p.Results.plot;
    band = p.Results.band;
    histon = p.Results.hist;
    topFreq = p.Results.topfreq;

    global alpha;
    alpha = p.Results.alpha;

    Suffix = '.mat';
    SuffixLen = size(Suffix, 2) - 1;

    [path, name, ext] = fileparts(basename);
    dataFolder = fullfile(path, name);
    d = dir(fullfile(path, name));
    loop = size(d, 1);

    normalPyrS = [];
    normalIntS = [];
    normalPyrnS = [];
    normalIntnS = [];
    pdPyrS = [];
    pdIntS = [];
    pdPyrnS = [];
    pdIntnS = [];

    possibleId = [];

    for i = 1:loop

        if length(d(i).name) > SuffixLen
            possibleId = [possibleId i];
        end

    end

    c = 1;

    for i = possibleId

        if strncmp(d(i).name( end - SuffixLen:end), Suffix, SuffixLen)
        filename = fullfile(dataFolder, d(i).name);
        fprintf('loading %s\n', filename);
        load(filename, 'phaseHistPyr', 'phaseHistInt', 'PyrIntList', 'PyrIntListStim');

        switch lower(stim)
            case 'both',
                normalPyrS = cat(3, normalPyrS, phaseHistPyr(:, :, PyrIntList{1}));
                normalIntS = cat(3, normalIntS, phaseHistInt(:, :, PyrIntList{2}));
                pdPyrS = cat(3, pdPyrS, phaseHistPyr(:, :, PyrIntList{3}));
                pdIntS = cat(3, pdIntS, phaseHistInt(:, :, PyrIntList{4}));
                TcharNormal = ['no stim'];
                TcharPD = ['no stim'];

            case 'pd',

                spdPyrS = phaseHistPyr(:, :, PyrIntListStim{3});
                spdIntS = phaseHistInt(:, :, PyrIntListStim{4});
                pdPyrS = cat(3, pdPyrS, phaseHistPyr(:, :, PyrIntListStim{3}));
                pdIntS = cat(3, pdIntS, phaseHistInt(:, :, PyrIntListStim{4}));

                nstim = setdiff(PyrIntList{3}, PyrIntListStim{3});
                spdPyrnS = phaseHistPyr(:, :, nstim);
                pdPyrnS = cat(3, pdPyrnS, phaseHistPyr(:, :, nstim));
                nstim = setdiff(PyrIntList{4}, PyrIntListStim{4});
                spdIntnS = phaseHistInt(:, :, nstim);
                pdIntnS = cat(3, pdIntnS, phaseHistInt(:, :, nstim));

                TcharNormal = ['PD stim'];
                TcharPD = ['PD stim'];
                RcharNormal = ['PD rec'];
                RcharPD = ['normal rec'];

            case 'normal',
                snormalPyrS = phaseHistPyr(:, :, PyrIntListStim{1});
                snormalIntS = phaseHistInt(:, :, PyrIntListStim{2});
                normalPyrS = cat(3, normalPyrS, phaseHistPyr(:, :, PyrIntListStim{1}));
                normalIntS = cat(3, normalIntS, phaseHistInt(:, :, PyrIntListStim{2}));

                nstim = setdiff(PyrIntList{1}, PyrIntListStim{1});
                snormalPyrnS = phaseHistPyr(:, :, nstim);
                normalPyrnS = cat(3, normalPyrnS, phaseHistPyr(:, :, nstim));
                nstim = setdiff(PyrIntList{2}, PyrIntListStim{2});
                snormalIntnS = phaseHistInt(:, :, nstim);
                normalIntnS = cat(3, normalIntnS, phaseHistInt(:, :, nstim));

                TcharNormal = ['normal stim'];
                TcharPD = ['normal stim'];
                RcharNormal = ['normal rec'];
                RcharPD = ['PD rec'];

            case 'stim',
                normalPyrS = cat(3, normalPyrS, phaseHistPyr(:, :, PyrIntListStim{1}));
                normalIntS = cat(3, normalIntS, phaseHistInt(:, :, PyrIntListStim{2}));
                pdPyrS = cat(3, pdPyrS, phaseHistPyr(:, :, PyrIntListStim{3}));
                pdIntS = cat(3, pdIntS, phaseHistInt(:, :, PyrIntListStim{4}));

                TcharNormal = ['norrmal stim'];
                TcharPD = ['PD stim'];
                RcharNormal = ['normal rec'];
                RcharPD = ['PD rec'];

            case 'nonstim',
                nstim = setdiff(PyrIntList{1}, PyrIntListStim{1});
                normalPyrS = cat(3, normalPyrS, phaseHistPyr(:, :, nstim));
                nstim = setdiff(PyrIntList{2}, PyrIntListStim{2});
                normalIntS = cat(3, normalIntS, phaseHistInt(:, :, nstim));
                nstim = setdiff(PyrIntList{3}, PyrIntListStim{3});
                pdPyrS = cat(3, pdPyrS, phaseHistPyr(:, :, nstim));
                nstim = setdiff(PyrIntList{4}, PyrIntListStim{4});
                pdIntS = cat(3, pdIntS, phaseHistInt(:, :, nstim));

                TcharNormal = ['normal stim'];
                TcharPD = ['PD stim'];
                RcharNormal = ['PD rec'];
                RcharPD = ['normal rec'];

        end

        if ~strcmp(plt, 'all')

            if isempty(strfind(stim, 'stim'))

                figure;

                switch lower(stim)
                    case 'pd',

                        subplot(2, 2, 1);
                        plotGainMap(spdPyrS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                        title(['PYR:' TcharNormal '- ' RcharNormal]);

                        subplot(2, 2, 2);
                        plotGainMap(spdIntS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                        title(['INT:' TcharNormal '- ' RcharNormal]);

                        subplot(2, 2, 3);
                        plotGainMap(spdPyrnS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                        title(['PYR:' TcharPD '- ' RcharPD]);

                        subplot(2, 2, 4);
                        plotGainMap(spdIntnS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                        title(['INT:' TcharPD '- ' RcharPD]);

                    case 'normal',
                        subplot(2, 2, 1);
                        plotGainMap(snormalPyrS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                        title(['PYR:' TcharNormal '- ' RcharNormal]);

                        subplot(2, 2, 2);
                        plotGainMap(snormalIntS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                        title(['INT:' TcharNormal '- ' RcharNormal]);

                        subplot(2, 2, 3);
                        plotGainMap(snormalPyrnS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                        title(['PYR:' TcharPD '- ' RcharPD]);

                        subplot(2, 2, 4);
                        plotGainMap(snormalIntnS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                        title(['INT:' TcharPD '- ' RcharPD]);
                end

                title(filename);
            end

        end

    end

end

if strcmp(plt, 'all')
    %pd or normal

    if isempty(strfind(stim, 'stim'))

        switch lower(stim)
            case 'pd',

                if ~strcmp(disp, 'both')
                    figure;
                    subplot(2, 2, 1);
                    mphpdPyrS = plotGainMap(pdPyrS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['PYR:' TcharNormal '- ' RcharNormal]);

                    subplot(2, 2, 2);
                    mphpdIntS = plotGainMap(pdIntS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['INT:' TcharNormal '- ' RcharNormal]);
                    X = mphpdPyrS;
                    Y = mphpdIntS;

                    subplot(2, 2, 3);
                    plotGainMap(pdPyrnS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['PYR:' TcharPD '- ' RcharPD]);

                    subplot(2, 2, 4);
                    plotGainMap(pdIntnS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['INT:' TcharPD '- ' RcharPD]);

                else
                    figure;
                    plotGainMap(pdPyrS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['PYR:' TcharNormal '- ' RcharNormal]);

                    figure;
                    plotGainMap(pdIntS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['INT:' TcharNormal '- ' RcharNormal]);

                    figure;
                    plotGainMap(pdPyrnS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['PYR:' TcharPD '- ' RcharPD]);

                    figure;
                    plotGainMap(pdIntnS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['INT:' TcharPD '- ' RcharPD]);

                end

            case 'normal',

                if ~strcmp(disp, 'both')
                    figure;
                    subplot(2, 2, 1);
                    mphnormalPyrS = plotGainMap(normalPyrS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['PYR:' TcharNormal '- ' RcharNormal]);

                    subplot(2, 2, 2);
                    mphnormalIntS = plotGainMap(normalIntS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['INT:' TcharNormal '- ' RcharNormal]);

                    X = mphnormalPyrS;
                    Y = mphnormalIntS;

                    subplot(2, 2, 3);
                    plotGainMap(normalPyrnS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['PYR:' TcharPD '- ' RcharPD]);

                    subplot(2, 2, 4);
                    plotGainMap(normalIntnS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['INT:' TcharPD '- ' RcharPD]);

                else

                    figure;
                    plotGainMap(normalPyrS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['PYR:' TcharNormal '- ' RcharNormal]);

                    figure;
                    plotGainMap(normalIntS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['INT:' TcharNormal '- ' RcharNormal]);

                    figure;
                    plotGainMap(normalPyrnS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['PYR:' TcharPD '- ' RcharPD]);

                    figure;
                    plotGainMap(normalIntnS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
                    title(['INT:' TcharPD '- ' RcharPD]);

                end

        end

        %statistical test for Gain to phase/freq relationship
        if (strcmp(band, 'freq') | strcmp(band, 'phase'))
            value = [X Y];
            group = [repmat((1:size(X, 1))', 1, size(X, 2)) repmat((1:size(Y, 1))', 1, size(Y, 2))];
            %group = [ones(size(X)) ones(size(Y))*2];
            %save test.mat X Y
            for q = 1:size(X, 1)
                [p, tbl, stats] = kruskalwallis([X(q, :) Y(q, :)], [zeros(1, size(X, 2)) ones(1, size(Y, 2))], 'off');

                if p < 0.01 / size(X, 1)
                    q
                    p
                end

            end

            phaseC = [-pi:pi / 4:pi];
            phaseC = phaseC(2:end);
            value = value';
            t1 = array2table(value);
            neuron = [zeros(1, size(X, 2)) ones(1, size(Y, 2))]';
            neuron = nominal(neuron);
            t2 = table(value(:, 1), neuron, 'variableNames', {'value1', 'neuron'});

            t = outerjoin(t1, t2, 'mergekeys', 1);
            within = [1:10]';

            rm = fitrm(t, 'value1-value10~neuron', 'withindesign', within);
            rtbl = ranova(rm, 'withinmodel', 'Time');
            
            if rtbl.pValue(5) < alpha
                c = multcompare(rm, 'neuron', 'by', 'Time');
                phaseC(c.Time(c.pValue(1:2:20) < alpha))
            end

        end

    else

        figure;
        subplot(2, 2, 1);
        plotGainMap(normalPyrS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
        title(['PYR:' TcharNormal '- ' RcharNormal]);

        subplot(2, 2, 2);
        plotGainMap(normalIntS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
        title(['INT:' TcharNormal '- ' RcharNormal]);

        subplot(2, 2, 3);
        plotGainMap(pdPyrS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
        title(['PYR:' TcharPD '- ' RcharPD]);

        subplot(2, 2, 4);
        plotGainMap(pdIntS, 'display', disp, 'band', band, 'hist', histon, 'topFreq', topFreq);
        title(['INT:' TcharPD '- ' RcharPD]);

    end

    colormap jet;
end
