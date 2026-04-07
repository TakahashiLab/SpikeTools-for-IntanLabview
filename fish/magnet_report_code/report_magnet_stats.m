function OUT = report_magnet_stats(csvPath, conditions, winStartSec, winLenSec, varargin)
% report_magnet_stats  v5.4 (circular pairwise 可視化 修正版)
% - PNGは既定で出力しない（PDFのみ）

%% Options
p = inputParser;
addParameter(p, 'outdir', 'magnet_report', @(x)ischar(x)||isstring(x));
addParameter(p, 'alpha', 0.05, @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'adjust', 'holm', @(x) any(strcmpi(x,{'holm','fdr'})));
addParameter(p, 'winlist', [], @(x) (isnumeric(x)&&isempty(x)) || (isnumeric(x)&&size(x,2)==2));
addParameter(p, 'doBootstrapDeltaR', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'bootstrapB', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p, 'useParallel', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'pdfRasterRes', 300, @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'drawMeanVector', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'drawTargetOnRose', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'roseBinN', 30, @(x)isnumeric(x)&&isscalar(x)&&x>=4);
addParameter(p, 'roseWritePNG', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'perFishWritePNG', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'fishMeansWritePNG', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'doUniformityTest', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'doTargetTest', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'rotateControlDeg', 0, @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'rotateControlNames', "control", @(x)ischar(x)||isstring(x)||iscellstr(x));
addParameter(p, 'targetDirDeg', NaN, @(x)isnumeric(x));
addParameter(p, 'doPerFishPlots', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'perFishTiles', [3 4], @(x)isnumeric(x)&&numel(x)==2);
addParameter(p, 'annotatePerFishP', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'doPerFishRayleigh', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'doFishMeansPlots', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'fishMeansTiles', [2 3], @(x)isnumeric(x)&&numel(x)==2);
addParameter(p, 'doStatsRaincloud', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'rainTiles', [], @(x)isnumeric(x)&&(isempty(x)||(numel(x)==2)));
addParameter(p, 'rainWritePNG', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'doCircPairPlot', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'circPairTiles', [], @(x)isnumeric(x)&&(isempty(x)||(numel(x)==2)));
addParameter(p, 'circPairWritePNG', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'doLME', false, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'modelMethod', 'auto', @(x) any(strcmpi(string(x), {'auto','lme','lm','off'})));
addParameter(p, 'includeOrderEffect', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'includeOrderInteraction', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'includePrevConditionEffect', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'cohortFilter', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(p, 'cohortGroupMode', 'sequence', @(x) any(strcmpi(string(x), {'sequence','cutoff'})));
addParameter(p, 'includeCohortGroupEffect', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'includeCohortGroupInteraction', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'cohortSplitCutoff', 12, @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'cohortGroupNames', {'cohort1to12','cohort13plus'}, @(x) iscell(x) && ~isempty(x));
addParameter(p, 'analyzeCohortGroupsSeparately', false, @(x)islogical(x)&&isscalar(x));
parse(p, varargin{:});
opt = p.Results;

if ischar(conditions) || isstring(conditions), conditions = cellstr(conditions); end
assert(numel(conditions)>=2, 'conditions は2つ以上を指定してください。');

if ~isempty(opt.winlist)
    if ~isfolder(opt.outdir), mkdir(opt.outdir); end
    OUT = struct(); OUT.conditions = string(conditions);
    OUT.winlist = opt.winlist; OUT.results = cell(size(opt.winlist,1),1);
    for i = 1:size(opt.winlist,1)
        ws = opt.winlist(i,1); wl = opt.winlist(i,2);
        subOutDir = fullfile(opt.outdir, sprintf('win_%gs_%gs', ws, finite_end(ws, wl)));
        res = report_magnet_stats_core(csvPath, conditions, ws, wl, opt, subOutDir);
        OUT.results{i} = res;
    end
    return
end

OUT = report_magnet_stats_core(csvPath, conditions, winStartSec, winLenSec, opt, string(opt.outdir));
if opt.analyzeCohortGroupsSeparately
    baseOpt = opt;
    baseOpt.analyzeCohortGroupsSeparately = false;
    names = string(opt.cohortGroupNames);
    ranges = [1, opt.cohortSplitCutoff; opt.cohortSplitCutoff + 1, inf];
    OUT.byCohortGroup = struct();
    for gi = 1:2
        subOpt = baseOpt;
        subOpt.cohortFilter = ranges(gi,:);
        subDir = fullfile(string(opt.outdir), char(names(gi)));
        OUT.byCohortGroup.(matlab.lang.makeValidName(char(names(gi)))) = ...
            report_magnet_stats_core(csvPath, conditions, winStartSec, winLenSec, subOpt, subDir);
    end
end
end
