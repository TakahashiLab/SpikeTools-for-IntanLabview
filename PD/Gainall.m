%
function [normalPyrS, normalIntS, pdPyrS, pdIntS] = Gainall(basename, varargin)
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
    p.addParamValue('cq',[25 .05], @isvector);

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
    CQs = p.Results.cq;

    global alpha;
    alpha = p.Results.alpha;

    Suffix = '.mat';
    SuffixLen = size(Suffix, 2) - 1;

    IDTh=CQs(1);
    LratioTh=CQs(2);

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

    possibleId = [];
    cellclassfile='cellclass.mat';
    for i = 1:loop

        if length(d(i).name) > SuffixLen & ~strcmp(d(i).name,cellclassfile)
            possibleId = [possibleId i];
        end

    end

    c = 1;

    for i = possibleId

        if strncmp(d(i).name( end - SuffixLen:end), Suffix, SuffixLen)
        filename = fullfile(dataFolder, d(i).name);
        fprintf('loading %s\n', filename);
        load(filename, 'phaseHistPyr', 'phaseHistInt', 'PyrIntList', 'PyrIntListStim', 'phaseHistPyrCtrl', 'phaseHistIntCtrl','pyr','interneuron','cq');
            %phaseHistPyr=phaseHistPyrCtrl;
            %phaseHistInt=phaseHistIntCtrl;
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
                normalPyr2=intersect(normalPyr,PyrIntList{9});
                normalInt2=intersect(normalInt,PyrIntList{10});
                pdPyr2=intersect(pdPyr,PyrIntList{9});
                pdInt2=intersect(pdInt,PyrIntList{10});
            case 'pv',
                normalPyr2=[];
                normalInt2=intersect(normalInt,PyrIntList{11});
                pdPyr2=[];
                pdInt2=intersect(pdInt,PyrIntList{11});
            case 'ccpv',
                normalPyr2=intersect(normalPyr,PyrIntList{9});
                normalInt2=intersect(normalInt,union(PyrIntList{10},PyrIntList{11}));
                pdPyr2=intersect(pdPyr,PyrIntList{9});
                pdInt2=intersect(pdInt,union(PyrIntList{10},PyrIntList{11}));
            case 'all',
                normalPyr2=[];
                normalInt2=[];
                pdPyr2=[];
                pdInt2=[];
        end
    
         %cell classification
         cq(isnan(cq(:,1)),1)=0;
         PYR=find(cq(pyr,1)>=IDTh & cq(pyr,2)<LratioTh);
    
        normalPyr=intersect(PYR,normalPyr);
         pdPyr=intersect(PYR,pdPyr);
   
         INTERNEURON=find(cq(interneuron,1)>=IDTh & cq(interneuron,2)<LratioTh);
         normalInt=intersect(INTERNEURON,normalInt);
         pdInt=intersect(INTERNEURON,pdInt);

        normalPyr=union(normalPyr,normalPyr2);
        normalInt=union(normalInt,normalInt2);
        pdPyr=union(pdPyr,pdPyr2);
        pdInt=union(pdInt,pdInt2);

     

        %
        normalPyrS = cat(3, normalPyrS, phaseHistPyr(:, :, normalPyr));
        normalIntS = cat(3, normalIntS, phaseHistInt(:, :, normalInt));
        normalPyrSCtrl = cat(3, normalPyrSCtrl, phaseHistPyrCtrl(:, :, normalPyr));
        normalIntSCtrl = cat(3, normalIntSCtrl, phaseHistIntCtrl(:, :, normalInt));
        
        pdPyrS = cat(3, pdPyrS, phaseHistPyr(:, :, pdPyr));
        pdIntS = cat(3, pdIntS, phaseHistInt(:, :, pdInt));
        pdPyrSCtrl = cat(3, pdPyrSCtrl, phaseHistPyrCtrl(:, :, pdPyr));
        pdIntSCtrl = cat(3, pdIntSCtrl, phaseHistIntCtrl(:, :, pdInt));

        if ~strcmp(plt, 'all')

            if isempty(strfind(stim, 'stim'))

                figure;

                switch lower(stim)
                    case 'pd',

                        subplot(2, 2, 1);
                        plotGainMap(pdPyrS, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
                        title(['PYR:' TcharNormal '- ' RcharNormal]);

                        subplot(2, 2, 2);
                        plotGainMap(pdIntS, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
                        title(['INT:' TcharNormal '- ' RcharNormal]);

                        subplot(2, 2, 3);
                        plotGainMap(pdPyrnS, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
                        title(['PYR:' TcharPD '- ' RcharPD]);

                        subplot(2, 2, 4);
                        plotGainMap(pdIntnS, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
                        title(['INT:' TcharPD '- ' RcharPD]);

                    case 'normal',
                        subplot(2, 2, 1);
                        plotGainMap(snormalPyrS, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
                        title(['PYR:' TcharNormal '- ' RcharNormal]);

                        subplot(2, 2, 2);
                        plotGainMap(snormalIntS, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
                        title(['INT:' TcharNormal '- ' RcharNormal]);

                        subplot(2, 2, 3);
                        plotGainMap(snormalPyrnS, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
                        title(['PYR:' TcharPD '- ' RcharPD]);

                        subplot(2, 2, 4);
                        plotGainMap(snormalIntnS, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
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

        if 0
        %if strcmp(disp, 'both')
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
            mphpdPyrS=plotGainMap(pdPyrS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['PYR:' TcharPD '- ' RcharPD]);

            figure;
            Ctrl{1} = pdIntSCtrl;
            mphpdIntS=plotGainMap(pdIntS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['INT:' TcharPD '- ' RcharPD]);

            p=mult_comp_perm_corr(mphnormalPyrS,mphpdPyrS);
            
        else
           
            figure;
            subplot(2, 2, 1);
            Ctrl{1} = normalPyrSCtrl;
            mphnormalPyrS = plotGainMap(normalPyrS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['PYR:' TcharNormal '- ' RcharNormal]);

            subplot(2, 2, 2);
            Ctrl{1} = normalIntSCtrl;
            mphnormalIntS = plotGainMap(normalIntS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['INT:' TcharNormal '- ' RcharNormal]);
            X = mphnormalPyrS;
            Y = mphnormalIntS;

            subplot(2, 2, 3);
            Ctrl{1} = pdPyrSCtrl;
            mphpdPyrS=plotGainMap(pdPyrS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['PYR:' TcharPD '- ' RcharPD]);

            subplot(2, 2, 4);
            Ctrl{1} = pdIntSCtrl;
            mphpdIntS=plotGainMap(pdIntS, 'control', Ctrl, 'display', disp, 'band', band, 'xyaxis', xyaxis, 'hist', histon, 'topFreq', topFreq);
            title(['INT:' TcharPD '- ' RcharPD]);
            
         
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
