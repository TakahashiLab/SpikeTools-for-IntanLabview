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
    cellclassfile='cellclass.mat';

    [path, name, ext] = fileparts(basename);
    dataFolder = fullfile(path, name);
    d = dir(fullfile(path, name));
    loop = size(d, 1);

    possibleId = [];

    for i = 1:loop

        if length(d(i).name) > SuffixLen & ~strcmp(d(i).name,cellclassfile)
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
GU = find(CQs(:, 2) < .1);

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
idx(~(p(:, 1) < 0.05 | p(:, 2) < 0.05)) = 3;

%check the cell type consistency 
idxGU=idx(find(GU(ind)));
intN=1;
pyrN=2;
if sum(idxGU(PIs(ind)==2)==2) < sum(idxGU(PIs(ind)==2)==1)
    intN=2;
    pyrN=1;
end

gscatter(Xnew(:, 1), Xnew(:, 2), idx);

save(cellclassfile,'GMModel','X','Xnew','intN','pyrN');
for i = possibleId
    if strncmp(d(i).name( end - SuffixLen:end), Suffix, SuffixLen)
        filename = fullfile(dataFolder, d(i).name);
        fprintf('loading %s\n', filename);
        load(filename, 'tp', 'sw', 'pi', 'cq');
        X = [tp sw];
        [idx, ~, p] = cluster(GMModel, X);
        idx(~(p(:, 1) < 0.05 | p(:, 2) < 0.05)) = 3;
        [~,name]=fileparts(d(i).name);
        execbuf=[name '= idx;'];
        eval(execbuf);
       
        save(cellclassfile,name,'-append');
      
   end
end

return
