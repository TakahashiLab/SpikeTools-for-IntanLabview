function mdl = fitlm_safe(T, formula)
y = strtrim(formula(1:find(formula=='~',1)-1));
good = true(height(T),1);
if ismember(y, T.Properties.VariableNames)
    if isnumeric(T.(y))
        good = good & isfinite(T.(y));
    else
        good = good & ~ismissing(T.(y));
    end
end
vars = {'Condition','Cohort','Fish','OrderInSession','PrevCondition','CohortGroup'};
for k = 1:numel(vars)
    if ismember(vars{k}, T.Properties.VariableNames)
        if isnumeric(T.(vars{k}))
            good = good & isfinite(T.(vars{k}));
        else
            good = good & ~ismissing(T.(vars{k}));
        end
    end
end
TT = T(good,:);
if isempty(TT)
    mdl = [];
    return;
end
try
    mdl = fitlm(TT, formula);
catch ME
    warning('fitlm_safe:FitFailed', 'Skipping LM for formula "%s": %s', formula, ME.message);
    mdl = [];
end
end
