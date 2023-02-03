%dn: 2018,2019,2020
function [phaseHistPyr, phaseHistInt, PyrIntList, PyrIntListStim, FRs, TPs, SWs, PIs, phaseHistPyrCtrl, phaseHistIntCtrl, CQs] = contPD(dn, an, varargin)

    p = inputParser;
    p.addParamValue('method', 'gainmap', @ischar);
    p.addParamValue('wavetype', 'chirp', @ischar);
    p.addParamValue('intensity', 5, @isnumeric);
    p.addParamValue('cellclass', 0, @isnumeric); %cell classification mode
    p.addParamValue('localcell', 0, @isnumeric); %cells only monitored from an optogenetic tetrode

    p.parse(varargin{:});
    method = p.Results.method;
    waveType = p.Results.wavetype;
    Intensity = p.Results.intensity;
    cellClass = p.Results.cellclass;
    localCell = p.Results.localcell;

    alpha = 0.05;
    [homePath, dataPath] = PDdataC(dn, an);
    TimingData = 'timing.mat';
    EnsembleData = 'ensemble.mat';
    ClassifyData = 'rensemble.mat';
    details = dataPath;

    loop = size(dataPath, 1);
    params = cell(loop, 2);

    phaseHistPyr = [];
    phaseHistInt = [];
    phaseHistPyrCtrl = [];
    phaseHistIntCtrl = [];
    offSetPyr = 0;
    offSetInt = 0;

    normalPyr = [];
    normalInt = [];
    PDpyr = [];
    PDint = [];
    ledPyr = [];
    ledInt = [];
    tagPyr = [];
    tagInt = [];

    CCPyr = [];
    CCInt = [];
    PVInt = [];

    normalPyrStim = [];
    normalIntStim = [];
    PDpyrStim = [];
    PDintStim = [];

    PyrIntList = [];
    PyrIntListStim = [];

    FRs = [];
    TPs = [];
    SWs = [];
    PIs = [];
    CQs = [];

    if strcmp(method, 'gainmap') | strcmp(method, 'coherence')
        method2 = method;
        method = 'gainOrCohere';
    end

    if strcmp(method, 'gainOrCohere')
        filename = an;
        load('cellclass.mat', filename, 'intN', 'pyrN');
        eval(['pi=' filename ';']);
        cPI = 1;
    end

    for i = 1:loop

        if dataPath{i, 7} == Intensity
            fprintf('analyzing %s\n', dataPath{i, 1});
            [~, name] = fileparts(dataPath{i, 1});

            switch lower(method)
                case 'behavior',
                    loadname = fullfile(homePath, dataPath{i, 1}, TimingData);
                    load(loadname, 'segPara', 'TrialT');

                    if dataPath{i, 4} <= 4 %LFP reference was on the left hemisphere

                        if dataPath{i, 4} == dataPath{i, 3}
                            LEDout = [LEDout; cellfun(@median, segPara(:, 2))']; %1:speed, 2:maxdist
                        else
                            LFPoutL = [LFPoutL; cellfun(@median, segPara(:, 2))']; %1:speed, 2:maxdist
                        end

                    else %Optical fiber was implanted in the right hemisphere

                        if dataPath{i, 4} == dataPath{i, 3}
                            LEDout = [LEDout; cellfun(@median, segPara(:, 2))']; %1:speed, 2:maxdist
                        else
                            LFPoutR = [LFPoutR; cellfun(@median, segPara(:, 2))']; %1:speed, 2:maxdist
                        end

                    end

                case 'cellclassify',
                    loadname = fullfile(homePath, dataPath{i, 1}, ClassifyData);

                    if exist(loadname, 'file')
                        load(loadname, 'rensemble');
                        TPs = [TPs; cell2mat(rensemble(:, 4))];
                        SWs = [SWs; cell2mat(rensemble(:, 5))];

                        pi = zeros(size(rensemble, 1), 1);
                        %interneuron (corr)

                        pi(dataPath{i, 8}) = 1;
                        %pyramidal cell(corr)
                        pi(dataPath{i, 9}) = 2;
                        %PV (tag)
                        pi(dataPath{i, 10}) = 3;

                        PIs = [PIs; pi];

                        loadname = fullfile(homePath, dataPath{i, 1}, EnsembleData);
                        load(loadname, 'IS', 'LR');
                        CQs = [CQs; IS LR];
                    end

                case 'gainorcohere',

                    loadname = fullfile(homePath, dataPath{i, 1}, TimingData);
                    load(loadname, 'preSilentTime', 'postSilentTime', 'chirpSeqTime', 'noiseSeqTime');
                    loadname = fullfile(homePath, dataPath{i, 1}, EnsembleData);
                    load(loadname, 'ensemble', 'an', 'en');
                    
                    piCP = pi(cPI:cPI + size(ensemble, 1) - 1);
                    cPI = cPI + size(ensemble, 1);
                    fprintf('real interneuron\n');

                    interneuron = union(dataPath{i, 8}, dataPath{i, 10}); % % % % %CC+tag
                    interneuronPV = dataPath{i, 10};
                    interneuronCC = dataPath{i, 8};
                    interneuron = union(interneuron, find(piCP == intN));

                    pyrTag = setdiff(dataPath{i, 11}, interneuron); % %conflict
                    pyrCC = union(dataPath{i, 9}, pyrTag); % % %CC+tag

                    pyr = union(pyrCC, find(piCP == pyrN));

                    if dataPath{i, 3} <= 4
                        lcq = 1:4;
                    else
                        lcq = 5:8;
                    end

                    lc0 = find(tetrodeMap(an, en) == dataPath{i, 3}); % units detected from a tetrode with an optical fiber
                    lcR = []; %units detected from a tetrode, whose spikes tagged with optical signals

                    for q = lcq
                        cand = find(tetrodeMap(an, en) == q);

                        if ~isempty(intersect(cand, dataPath{i, 10}))
                            lcR = [lcR cand];
                        end

                    end

                    %interneuron = intersect(interneuron, lc);
                    %pyr = intersect(pyr, lc);

                    switch lower(waveType)
                        case 'chirp',
                            SeqTime = chirpSeqTime;
                        case 'noise',
                            SeqTime = noiseSeqTime;
                    end

                    if cellClass
                        %classify pyr and int
                        [pyr, interneuron, fr, tp, sw] = checkFR(ensemble);
                        FRs = [FRs; fr];
                        TPs = [TPs; tp];
                        SWs = [SWs; sw];
                        PIs = [PIs; sparse(dataPath{i, 8}, 1, 1, length(ensemble), 1)];
                    else

                        if ~isempty(pyr) | ~isempty(interneuron)
                            loadname = fullfile(homePath, dataPath{i, 1}, [name 'Event.mat']);
                            load(loadname, 'event');

                            if strcmp(method2, 'gainmap')
                                %calc gain map
                                [pPyr, pInt, pPyrCtrl, pIntCtrl] = batchGainMap(event, ensemble(:, 3), SeqTime, 1, [preSilentTime; postSilentTime], pyr, interneuron);
                            elseif strcmp(method2, 'coherence')
                                %[Cs,phis,fs]=batchMTS(LFP,ensemble,chirpSeqTime,tetNum);
                                %coherence LED and neuronal ensemble
                                [Cs, phis, fs, confCs] = batchMTS(event, ensemble(:, 3), SeqTime, 1);
                                %plotMTS(Cs,phis,fs,pyr,interneuron,confCs);
                            end

                            %cell classification
                            normal = find(tetrodeMap(an, en) <= 4);

                            %putative pyramidal cell
                            if ~isempty(pyr)
                                %pyramidal

                                [~, ind] = intersect(pyr, 1:normal(end));
                                normalPyr = [normalPyr; ind + offSetPyr];

                                if dataPath{i, 3} <= 4
                                    normalPyrStim = [normalPyrStim; ind + offSetPyr];
                                end

                                [~, ind] = intersect(pyr, pyrCC);
                                CCPyr = [CCPyr; ind + offSetPyr];

                                if ~isempty(lc0)
                                    [~, ind] = intersect(pyr, lc0);
                                    ledPyr = [ledPyr; ind + offSetPyr];
                                end

                                if ~isempty(lcR)
                                    [~, ind] = intersect(pyr, lcR);
                                    tagPyr = [tagPyr; ind + offSetPyr];
                                end

                                [~, ind] = setdiff(pyr, 1:normal(end));
                                PDpyr = [PDpyr; ind + offSetPyr];

                                if dataPath{i, 3} >= 5
                                    PDpyrStim = [PDpyrStim; ind + offSetPyr];
                                end

                            end

                            %interneuron or PV
                            if ~isempty(interneuron)
                                %interneuron
                                [~, ind] = intersect(interneuron, 1:normal(end));
                                normalInt = [normalInt; ind + offSetInt];

                                if dataPath{i, 3} <= 4
                                    normalIntStim = [normalIntStim; ind + offSetInt];
                                end

                                [~, ind] = intersect(interneuron, interneuronCC);
                                CCInt = [CCInt; ind + offSetInt];

                                [~, ind] = intersect(interneuron, interneuronPV);
                                PVInt = [PVInt; ind + offSetInt];

                                if ~isempty(lc0)
                                    [~, ind] = intersect(interneuron, lc0);
                                    ledInt = [ledInt; ind + offSetInt];
                                end

                                if ~isempty(lcR)
                                    [~, ind] = intersect(interneuron, lcR);
                                    tagInt = [tagInt; ind + offSetInt];
                                end

                                [~, ind] = setdiff(interneuron, 1:normal(end));
                                PDint = [PDint; ind + offSetInt];

                                if dataPath{i, 3} >= 5
                                    PDintStim = [PDintStim; ind + offSetInt];
                                end

                            end

                            offSetPyr = offSetPyr + length(pyr);
                            offSetInt = offSetInt + length(interneuron);

                            if strcmp(method2, 'gainmap')
                                phaseHistPyr = cat(3, phaseHistPyr, pPyr);
                                phaseHistInt = cat(3, phaseHistInt, pInt);

                                phaseHistPyrCtrl = cat(3, phaseHistPyrCtrl, pPyrCtrl);
                                phaseHistIntCtrl = cat(3, phaseHistIntCtrl, pIntCtrl);
                            elseif strcmp(method2, 'coherence')
                                phaseHistPyr = cat(1, phaseHistPyr, Cs);
                                phaseHistInt = cat(1, phaseHistInt, phis);

                                phaseHistPyrCtrl = cat(1, phaseHistPyrCtrl, fs);
                                phaseHistIntCtrl = cat(1, phaseHistIntCtrl, confCs);

                            end

                        end

                    end

            end

        end

    end

    if ~cellClass
        PyrIntList{1} = normalPyr;
        PyrIntList{2} = normalInt;
        PyrIntList{3} = PDpyr;
        PyrIntList{4} = PDint;
        PyrIntList{5} = ledPyr;
        PyrIntList{6} = ledInt;
        PyrIntList{7} = tagPyr;
        PyrIntList{8} = tagInt;
        PyrIntList{9} = CCPyr;
        PyrIntList{10} = CCInt;
        PyrIntList{11} = PVInt;

        PyrIntListStim{1} = normalPyrStim;
        PyrIntListStim{2} = normalIntStim;
        PyrIntListStim{3} = PDpyrStim;
        PyrIntListStim{4} = PDintStim;
    end

    return;
