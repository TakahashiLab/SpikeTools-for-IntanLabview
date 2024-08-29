%
function [TPs, SWs, PIs] = Cellall(basename, varargin)
    p = inputParser;
    p.addParamValue('proc', 'initial', @ischar);
    p.addParamValue('alpha', 0.05, @isscalar);
    p.parse(varargin{:});

    proc = p.Results.proc;

    global alpha;
    alpha = p.Results.alpha;

    Suffix = '.mat';
    SuffixLen = size(Suffix, 2) - 1;
    cellclassfile = 'cellclass.mat';
    LratioTh = 0.05; %0.05, 0.1
    IDTh = 15;
    pTh = 0.05; %0.05, 0.1

    [path, name, ext] = fileparts(basename);
    dataFolder = fullfile(path, name);
    d = dir(fullfile(path, name));
    loop = size(d, 1);

    possibleId = [];

    for i = 1:loop

        if length(d(i).name) > SuffixLen & ~strcmp(d(i).name, cellclassfile)
            possibleId = [possibleId i];
        end

    end

    if strcmp(proc, 'initial')
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
    GU = find(CQs(:, 2) < LratioTh & CQs(:, 1) > IDTh); %L-ratio, Isolation Distance

    TPs = TPs(GU);
    SWs = SWs(GU);
    PIs = PIs(GU);

    X = [TPs SWs];
    ind = find(PIs == 1 | PIs == 2);
    Xnew = X;
    X = X(ind, :);

    GMModel = fitgmdist(X, 2, 'sharedcovariance', true, 'covariancetype', 'full', 'start', 'randSample');
    y = PIs(ind) - 1;

    subplot(1, 2, 1)
    h = gscatter(X(:, 1), X(:, 2), y);
    hold on
    gmPDF = @(x, y) arrayfun(@(x0, y0) pdf(GMModel, [x0 y0]), x, y);

    g = gca;
    fcontour(gmPDF, [0 1.0 0 1.6])
    title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
    legend(h, 'Model 0', 'Model1')

    subplot(1, 2, 2);
    [idx, ~, p] = cluster(GMModel, Xnew);

    if ~(pTh < 0)
        idx(~(p(:, 1) < pTh | p(:, 2) < pTh)) = 3;
    end

    %check the cell type consistency

    idxGU = idx(find(GU(ind)));
    intN = 1;
    pyrN = 2;

    if sum(idxGU(PIs(ind) == 1) == 1) > sum(idxGU(PIs(ind) == 2) == 1)
        intN = 2;
        pyrN = 1;
    end

    gscatter(Xnew(:, 1), Xnew(:, 2), idx);

    save(cellclassfile, 'GMModel', 'X', 'Xnew', 'intN', 'pyrN');
else
    load(cellclassfile, 'GMModel');
end

%summary statistics
fprintf('Passed units=%d\n', length(GU));
fprintf('pyramidalc cell=%d, interneuron=%d\n', sum(idx == pyrN), sum(idx == intN));
fprintf('Total cell number=%d\n', sum(idx == pyrN | idx == intN));

for i = possibleId
    if strncmp(d(i).name( end - SuffixLen:end), Suffix, SuffixLen)
    filename = fullfile(dataFolder, d(i).name);
    fprintf('loading %s\n', filename);
    load(filename, 'tp', 'sw', 'pi', 'cq');
    X = [tp sw];
    [idx, ~, p] = cluster(GMModel, X);
    idx(~sum(p < 0.05, 2)) = 3;
    %BU = find(cq(:, 2) >= LratioTh | cq(:, 1) <= IDTh); %L-ratio, Isolation Distance
    %idx(BU) = 4;
    [~, name] = fileparts(d(i).name);
    execbuf = [name '= idx;'];
    eval(execbuf);

    save(cellclassfile, name, '-append');

end

end
return;
