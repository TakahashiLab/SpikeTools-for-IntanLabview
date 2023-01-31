%
function [TPs, SWs, PIs] = Cellall(basename, varargin)
    p = inputParser;
    p.addParamValue('proc', 'individual', @ischar);
    p.addParamValue('stimulation', 'both', @ischar);
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

    possibleId = [];

    for i = 1:loop

        if length(d(i).name) > SuffixLen
            possibleId = [possibleId i];
        end

    end

    c = 1;

    TPs = [];
    SWs = [];
    PIs = [];
    CQs = [];

    for i = possibleId
        if strncmp(d(i).name( end - SuffixLen:end), Suffix, SuffixLen)
        filename = fullfile(dataFolder, d(i).name);
        fprintf('loading %s\n', filename);
        load(filename, 'tp', 'sw', 'pi', 'cq');
        TPs = [TPs; tp];
        SWs = [SWs; sw];
        PIs = [PIs; pi];
        CQs = [CQs; cq];
    end

end

%good units
GU = find(CQs(:, 2) < 0.05);
TPs = TPs(GU);
SWs = SWs(GU);
PIs = PIs(GU);

if 1
    %plot(TPs, SWs, 'k.');
    hold on;

    plot(TPs(PIs == 1), SWs(PIs == 1), 'bo');
    plot(TPs(PIs == 2), SWs(PIs == 2), 'go');
    scatter(TPs(PIs == 3), SWs(PIs == 3), 'o', 'markerfacecolor', 'b');
end

if 0
    idx = kmeans([TPs SWs], 2);
    figure;

    hold on;
    plot(TPs(idx == 1), SWs(idx == 1), 'bo');
    plot(TPs(idx == 2), SWs(idx == 2), 'go');
elseif 0
    X = [TPs SWs];
    ind = find(PIs == 1 | PIs == 2);
    X = sparse(X(ind, :));
    y = PIs(ind);
    Mdl = fitclinear(X, y);
    label = predict(Mdl, X);

    figure
    scatter(X(label == 1, 1), X(label == 1, 2), 'r.');
    hold on
    scatter(X(label == 2, 1), X(label == 2, 2), 'g.');
    scatter(X(y == 1, 1), X(y == 1, 2), 'ro');
    scatter(X(y == 2, 1), X(y == 2, 2), 'go');
    hold off
end

return
