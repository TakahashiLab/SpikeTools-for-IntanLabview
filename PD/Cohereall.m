%
function [normalPyrS, normalIntS, pdPyrS, pdIntS] = Cohereall(basename, varargin)
    p = inputParser;
    p.addParamValue('proc', 'individual', @ischar);
    p.addParamValue('stimulation', 'both', @ischar);
    p.addParamValue('celltype', 'all', @ischar);
    p.addParamValue('display', 'ascend', @ischar);
    p.addParamValue('plot', 'all', @ischar);
    p.addParamValue('band', 'non', @ischar);
    p.addParamValue('hist', 'off', @ischar);
    p.addParamValue('topfreq', 120, @isnumeric);
    p.addParamValue('alpha', 0.05, @isnumeric);
    p.addParamValue('xyaxis', 'off', @ischar);

    p.parse(varargin{:});

    proc = p.Results.proc;
    stim = p.Results.stimulation;
    celltype = p.Results.celltype;
    disp = p.Results.display;
    plt = p.Results.plot;
    band = p.Results.band;
    histon = p.Results.hist;
    topFreq = p.Results.topfreq;
    xyaxis = p.Results.xyaxis;

    global alpha;
    alpha = p.Results.alpha;

    Suffix = '.mat';
    SuffixLen = size(Suffix, 2) - 1;

    [path, name, ext] = fileparts(basename);
    dataFolder = fullfile(path, name);
    d = dir(fullfile(path, name));
    loop = size(d, 1);

    normalPyr = [];
    normalInt = [];
    pdPyr = [];
    pdInt = [];

    normalPyrS = [];
    normalIntS = [];
    normalPyrnS = [];
    normalIntnS = [];
    pdPyrS = [];
    pdIntS = [];
    pdPyrnS = [];
    pdIntnS = [];

    normalPyrSCtrl = [];
    normalIntSCtrl = [];
    normalPyrnSCtrl = [];
    normalIntnSCtrl = [];
    pdPyrSCtrl = [];
    pdIntSCtrl = [];
    pdPyrnSCtrl = [];
    pdIntnSCtrl = [];

    normalPyrs = [];
    normalInts = [];
    pdPyrs = [];
    pdInts = [];
    offSet = zeros(1, 2);

    possibleId = [];
    cellclassfile = 'cellclass.mat';

    for i = 1:loop

        if length(d(i).name) > SuffixLen & ~strcmp(d(i).name, cellclassfile)
            possibleId = [possibleId i];
        end

    end

    c = 1;
    div=0;
    switch (disp)
    case 'ascend',
        div=1;
    case 'descend',
        div=2;
    case 'both',
        div=3;
    end
    

    for i = possibleId

        if strncmp(d(i).name( end - SuffixLen:end), Suffix, SuffixLen)
        filename = fullfile(dataFolder, d(i).name);
        fprintf('loading %s\n', filename);
        load(filename, 'phaseHistPyr', 'phaseHistInt', 'PyrIntList', 'PyrIntListStim', 'phaseHistPyrCtrl', 'phaseHistIntCtrl', 'pyr', 'interneuron');
        [~, fn] = fileparts(filename);

        if 0
            %from cellclass.mat
            eval(['PI=' fn ';']);
            pyr = find(PI == pyrN); %classification
            pyr = unique([pyr; cell2mat(PyrIntList([1 3 9])')]); %cc, pv
            interneuron = find(PI == intN); %classification
            interneuron = unique([interneuron; cell2mat(PyrIntList([2 4 10 11])')]);
            PyrIntList{1}
            PyrIntList{2}
        end

        switch lower(stim)
            case 'both',
                normalPyr = PyrIntList{1};
                normalInt = PyrIntList{2};
                pdPyr = PyrIntList{3};
                pdInt = PyrIntList{4};

                TcharNormal = ['no stim'];
                TcharPD = ['no stim'];

            case 'pd',

                pdPyr = PyrIntListStim{3};
                pdInt = PyrIntListStim{4};

                nstim = setdiff(PyrIntList{3}, PyrIntListStim{3});
                pdPyrnS = cat(3, pdPyrnS, phaseHistPyr(:, :, nstim));
                pdPyrnSCtrl = cat(3, pdPyrnSCtrl, phaseHistPyrCtrl(:, :, nstim));
                nstim = setdiff(PyrIntList{4}, PyrIntListStim{4});
                pdIntnS = cat(3, pdIntnS, phaseHistInt(:, :, nstim));
                pdIntnSCtrl = cat(3, pdIntnSCtrl, phaseHistIntCtrl(:, :, nstim));

                TcharNormal = ['PD stim'];
                TcharPD = ['PD stim'];
                RcharNormal = ['PD rec'];
                RcharPD = ['normal rec'];

            case 'normal',
                normalPyr = PyrIntListStim{1};
                normalInt = PyrIntListStim{2};

                nstim = setdiff(PyrIntList{1}, PyrIntListStim{1});
                normalPyrnS = cat(3, normalPyrnS, phaseHistPyr(:, :, nstim));
                normalPyrnSCtrl = cat(3, normalPyrnSCtrl, phaseHistPyrCtrl(:, :, nstim));

                nstim = setdiff(PyrIntList{2}, PyrIntListStim{2});
                normalIntnS = cat(3, normalIntnS, phaseHistInt(:, :, nstim));
                normalIntnSCtrl = cat(3, normalIntnSCtrl, phaseHistIntCtrl(:, :, nstim));

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

            case 'ledstim',
                normalPyr = setdiff(PyrIntList{5}, PyrIntListStim{3});
                normalInt = setdiff(PyrIntList{6}, PyrIntListStim{4});
                TcharNormal = ['normal led neighbor'];
                RcharNormal = [];

                pdPyr = setdiff(PyrIntList{5}, PyrIntListStim{1});
                pdInt = setdiff(PyrIntList{6}, PyrIntListStim{2});
                TcharPD = ['pd led neighbor'];
                RcharPD = [];

            case 'tagstim',
                normalPyr = setdiff(PyrIntList{7}, PyrIntListStim{3});
                normalInt = setdiff(PyrIntList{8}, PyrIntListStim{4});
                TcharNormal = ['normal tagging units'];
                RcharNormal = [];

                pdPyr = setdiff(PyrIntList{7}, PyrIntListStim{1});
                pdInt = setdiff(PyrIntList{8}, PyrIntListStim{2});

                TcharPD = ['pd tagging neighbor'];
                RcharPD = [];

            case 'ledtagstim',
                normalPyr = setdiff(union(PyrIntList{5}, PyrIntList{7}), PyrIntListStim{3});
                normalInt = setdiff(union(PyrIntList{6}, PyrIntList{8}), PyrIntListStim{4});
                TcharNormal = ['normal led tag'];
                RcharNormal = [];

                pdPyr = setdiff(union(PyrIntList{5}, PyrIntList{7}), PyrIntListStim{1});
                pdInt = setdiff(union(PyrIntList{6}, PyrIntList{8}), PyrIntListStim{2});
                TcharPD = ['pd led tag'];
                RcharPD = [];

            case 'farstim',
                normalPyr = setdiff(PyrIntListStim{1}, union(PyrIntList{5}, PyrIntList{7}));
                normalInt = setdiff(PyrIntListStim{2}, union(PyrIntList{6}, PyrIntList{8}));
                TcharNormal = ['normal stim'];
                RcharNormal = ['normal(ipsi) far'];

                pdPyr = setdiff(PyrIntListStim{3}, union(PyrIntList{5}, PyrIntList{7}));
                pdInt = setdiff(PyrIntListStim{4}, union(PyrIntList{6}, PyrIntList{8}));
                TcharPD = ['pd stim'];
                RcharPD = ['pd(ipsi) far'];

            case 'contrastim',
                normalPyr = setdiff(PyrIntList{1}, PyrIntListStim{1});
                normalInt = setdiff(PyrIntList{2}, PyrIntListStim{2});
                TcharNormal = ['pd stim '];
                RcharNormal = ['normal contra'];

                pdPyr = setdiff(PyrIntList{3}, PyrIntListStim{3});
                pdInt = setdiff(PyrIntList{4}, PyrIntListStim{4});
                TcharPD = ['normal stim'];
                RcharPD = ['pd contra'];

        end

        switch (lower(celltype))
            case 'cc',
                normalPyr = intersect(normalPyr, PyrIntList{9});
                normalInt = intersect(normalInt, PyrIntList{10});
                pdPyr = intersect(pdPyr, PyrIntList{9});
                pdInt = intersect(pdInt, PyrIntList{10});
            case 'pv',
                normalInt = intersect(normalInt, PyrIntList{11});
                pdInt = intersect(pdInt, PyrIntList{11});
            case 'ccpv',
                normalPyr = intersect(normalPyr, PyrIntList{9});
                normalInt = intersect(normalInt, union(PyrIntList{10}, PyrIntList{11}));
                pdPyr = intersect(pdPyr, PyrIntList{9});
                pdInt = intersect(pdInt, union(PyrIntList{10}, PyrIntList{11}));
            case 'all',

        end

        normal = [pyr(normalPyr)' interneuron(normalInt)'];
        
        normalPyrs = [normalPyrs (1:length(normalPyr)) + offSet(1)];
       
        normalInts = [normalInts (length(normalPyr) + 1:length(normal)) + offSet(1)];
        
        if ~isempty(normal)
            normalPyrS = cat(1, normalPyrS, phaseHistPyr(normal));
            normalIntS = cat(1, normalIntS, phaseHistInt(normal));
            normalPyrSCtrl = cat(1, normalPyrSCtrl, phaseHistPyrCtrl(normal));
            normalIntSCtrl = cat(1, normalIntSCtrl, phaseHistIntCtrl(normal));
        end

        pd = [pyr(pdPyr)' interneuron(pdInt)'];

        pdPyrs = [pdPyrs (1:length(pdPyr)) + offSet(2)];
        pdInts = [pdInts (length(pdPyr) + 1:length(pd)) + offSet(2)];

        if ~isempty(pd)
            pdPyrS = cat(1, pdPyrS, phaseHistPyr(pd));
            pdIntS = cat(1, pdIntS, phaseHistInt(pd));
            pdPyrSCtrl = cat(1, pdPyrSCtrl, phaseHistPyrCtrl(pd));
            pdIntSCtrl = cat(1, pdIntSCtrl, phaseHistIntCtrl(pd));
        end

        offSet = offSet + [length(normal) length(pd)];
    end

end


if strcmp(plt, 'all')
    %pd or normal

    if isempty(strfind(stim, 'stim'))

        switch lower(stim)
            case 'pd',

            case 'normal',
                pdPyrS = normalPyrS;
                pdPyrSCtrl = normalPyrSCtrl;
                pdPyrnS = normalPyrnS;
                pdPyrnSCtrl = normalPyrnSCtrl;
                pdIntS = normalIntS;
                pdIntSCtrl = normalIntSCtrl;
                pdIntnS = normalIntnS;
                pdIntnSCtrl = normalIntnSCtrl;
        end

        if strcmp(disp, 'both')
            figure;
            Ctrl{1} = pdPyrSCtrl;
            mphpdPyrS = plotGainMap(pdPyrS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['PYR:' TcharNormal '- ' RcharNormal]);

            figure;
            Ctrl{1} = pdIntSCtrl;
            mphpdIntS = plotGainMap(pdIntS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['INT:' TcharNormal '- ' RcharNormal]);
            X = mphpdPyrS;
            Y = mphpdIntS;

            figure;
            Ctrl{1} = pdPyrnSCtrl;
            plotGainMap(pdPyrnS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['PYR:' TcharPD '- ' RcharPD]);

            figure;
            Ctrl{1} = pdIntnSCtrl;
            plotGainMap(pdIntnS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['INT:' TcharPD '- ' RcharPD]);
        else
            figure;
            subplot(2, 2, 1);
            Ctrl{1} = pdPyrSCtrl;
            mphpdPyrS = plotGainMap(pdPyrS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['PYR:' TcharNormal '- ' RcharNormal]);

            subplot(2, 2, 2);
            Ctrl{1} = pdIntSCtrl;
            mphpdIntS = plotGainMap(pdIntS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['INT:' TcharNormal '- ' RcharNormal]);
            X = mphpdPyrS;
            Y = mphpdIntS;

            subplot(2, 2, 3);
            Ctrl{1} = pdPyrnSCtrl;
            plotGainMap(pdPyrnS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['PYR:' TcharPD '- ' RcharPD]);

            subplot(2, 2, 4);
            Ctrl{1} = pdIntnSCtrl;
            plotGainMap(pdIntnS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['INT:' TcharPD '- ' RcharPD]);
        end

    elseif strcmp(stim, 'ledstim') | strcmp(stim, 'tagstim') | strcmp(stim, 'ledtagstim') | strcmp(stim, 'farstim') | strcmp(stim, 'contrastim')
        if 0
        %if strcmp(disp, 'both')
            figure;

            Ctrl{1} = normalPyrSCtrl;
            mphnormalPyrS = plotGainMap(normalPyrS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['PYR:' TcharNormal '- ' RcharNormal]);

            figure;
            Ctrl{1} = normalIntSCtrl;
            mphnormalIntS = plotGainMap(normalIntS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['INT:' TcharNormal '- ' RcharNormal]);
            X = mphPyrS;
            Y = mphIntS;

            figure;
            Ctrl{1} = pdPyrSCtrl;
            mphpdPyrS = plotGainMap(pdPyrS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['PYR:' TcharPD '- ' RcharPD]);

            figure;
            Ctrl{1} = pdIntSCtrl;
            mphpdIntS = plotGainMap(pdIntS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['INT:' TcharPD '- ' RcharPD]);

            p = mult_comp_perm_corr(mphnormalPyrS, mphpdPyrS);

        else
            
            Cs = normalPyrS;
            phis = normalIntS;
            fs = normalPyrSCtrl;
            confCs = normalIntSCtrl;
         
            fslabel = [2:2:topFreq];
         
            [thetaIdxNormal,peakCohNormal,peakFreqNormal] = plotMTS(Cs, phis, fs, normalPyrs, normalInts, confCs, fslabel, topFreq,div);

            Cs = pdPyrS;
            phis = pdIntS;
            fs = pdPyrSCtrl;
            confCs = pdIntSCtrl;

            figure;
            [thetaIdxPD,peakCohPD,peakFreqPD] = plotMTS(Cs, phis, fs, pdPyrs, pdInts, confCs, fslabel, topFreq,div);

            compThetaIdx(thetaIdxNormal, thetaIdxPD,'theta index');
            compThetaIdx(peakCohNormal,peakCohPD,'peak coherence');
            compThetaIdx(peakFreqNormal,peakFreqPD,'peak frequency (Hz)');

        end

    else %stim/nonstim

        figure;
        subplot(2, 2, 1);
        plotGainMap(normalPyrS, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
        title(['PYR:' TcharNormal '- ' RcharNormal]);

        subplot(2, 2, 2);
        plotGainMap(normalIntS, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
        title(['INT:' TcharNormal '- ' RcharNormal]);

        subplot(2, 2, 3);
        plotGainMap(pdPyrS, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
        title(['PYR:' TcharPD '- ' RcharPD]);

        subplot(2, 2, 4);
        plotGainMap(pdIntS, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
        title(['INT:' TcharPD '- ' RcharPD]);

    end

    colormap jet;

    %statistical test for Gain to phase/freq relationship
    if 0
        %if (strcmp(band,  'freq') | strcmp(band, 'phase'))

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
        rtbl.pValue(5)

        if rtbl.pValue(5) < alpha
            c = multcompare(rm, 'neuron', 'by', 'Time');

            if strcmp(band, 'phase')
                fprintf('phase=%f\n', phaseC(c.Time(c.pValue(1:2:20) < alpha)));
            else
                fprintf('freq=%f\n', c.Time(c.pValue(1:2:20) < alpha) * 2)
            end

        end

    end

end
end
%%%%%%%%%%%55
function compThetaIdx(thetaIdxNormal, thetaIdxPD,yLABEL)
    x = [cell2mat(thetaIdxNormal) cell2mat(thetaIdxPD)];
    grp = [ones(size(thetaIdxNormal{1})) ones(size(thetaIdxNormal{2})) * 2 ones(size(thetaIdxPD{1})) * 3 ones(size(thetaIdxPD{2})) * 4];
    [p, tbl, stat] = kruskalwallis(x, grp,'off');
   
    [c,m]=multcompare(stat,'CriticalValueType','bonferroni','Display','off');
   c
    [avg,se]=grpstats(x,grp,["mean","sem"]);
    figure;
    bar(avg);
    hold on;
    errorbar([1:4],avg,se,'.');
    ylabel(yLABEL);
    set(gca,'xtick',1:4,'xticklabel',{'n-Pyr','n-PV','pd-Pyr','pd-PV'})
    
end
