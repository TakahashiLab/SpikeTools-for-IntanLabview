function fileMap = draw_lme_anova_tables(lme_cont, anova_cos, anova_sin, outdir, ws, wl, pngRes)
if nargin < 7, pngRes = 300; end

pdfPath = fullfile(outdir, 'LME_anova_tables.pdf');
if exist(pdfPath,'file'), try, delete(pdfPath); end, end

pages = {};
metricNames = fieldnames(lme_cont);
for i = 1:numel(metricNames)
    name = metricNames{i};
    A = [];
    if isfield(lme_cont.(name), 'anova')
        A = lme_cont.(name).anova;
    end
    if ~isempty(A)
        pages(end+1,:) = {sprintf('LME ANOVA: %s', name), table_to_text(A)}; %#ok<AGROW>
    end
end
if ~isempty(anova_cos)
    pages(end+1,:) = {'LME ANOVA: cosMag', table_to_text(anova_cos)}; %#ok<AGROW>
end
if ~isempty(anova_sin)
    pages(end+1,:) = {'LME ANOVA: sinMag', table_to_text(anova_sin)}; %#ok<AGROW>
end

if isempty(pages)
    fileMap = struct('pdf', "", 'png', "");
    return;
end

for i = 1:size(pages,1)
    f = figure('Color','w','Units','pixels','Position',[100 100 1200 900], 'Visible','off');
    ax = axes('Parent', f, 'Position', [0 0 1 1], 'Visible', 'off'); %#ok<LAXES>
    titleStr = sprintf('%s | Window %gs-%gs', pages{i,1}, ws, finite_end(ws, wl));
    text(0.03, 0.97, titleStr, 'Units','normalized', ...
        'FontName','Consolas', 'FontSize', 14, 'FontWeight','bold', ...
        'VerticalAlignment','top', 'Interpreter','none');
    text(0.03, 0.92, pages{i,2}, 'Units','normalized', ...
        'FontName','Consolas', 'FontSize', 10, ...
        'VerticalAlignment','top', 'Interpreter','none');
    finalize_page(f, pdfPath, "", pngRes, i > 1, false);
    close(f);
end

fileMap = struct('pdf', pdfPath, 'png', "");
end

function txt = table_to_text(T)
try
    txt = strtrim(evalc('disp(T)'));
catch
    txt = '[Could not render ANOVA table]';
end
if isempty(txt)
    txt = '[Empty ANOVA table]';
end
end
