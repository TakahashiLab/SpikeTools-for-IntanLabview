function mdl = fitglme_safe(T, formula, varargin)
y = strtrim(formula(1:find(formula=='~',1)-1));
need = {'Condition','Cohort','Fish'};
good = true(height(T),1);
if ismember(y, T.Properties.VariableNames)
    if isnumeric(T.(y)), good = good & isfinite(T.(y));
    else,               good = good & ~ismissing(T.(y));
    end
end
for k=1:numel(need)
    if ismember(need{k}, T.Properties.VariableNames)
        if isnumeric(T.(need{k})); good = good & isfinite(T.(need{k}));
        else;                       good = good & ~ismissing(T.(need{k}));
        end
    end
end
TT = T(good,:);
if isempty(TT), mdl = []; return; end
try
    mdl = fitglme(TT, formula, varargin{:});
catch ME
    warning('fitglme_safe:FitFailed', 'Skipping LME for formula "%s": %s', formula, ME.message);
    mdl = [];
end
end
