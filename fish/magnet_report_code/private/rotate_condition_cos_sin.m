function S = rotate_condition_cos_sin(S, rotateDeg, condNames)
% rotate_condition_cos_sin  Rotate cos/sin components for specified conditions by rotateDeg (degrees)
% Usage:
%   S = rotate_condition_cos_sin(S, 15, "control");
%   S = rotate_condition_cos_sin(S, -30, {"control","controlSea"});
%
% Only S.cosMag and S.sinMag are rotated. Other columns are unchanged.

if isempty(S), return; end
if ischar(condNames) || isstring(condNames)
    condNames = cellstr(string(condNames));
elseif iscellstr(condNames)
    % ok
else
    error('condNames must be char, string, or cellstr');
end

mask = false(height(S),1);
for k = 1:numel(condNames)
    mask = mask | strcmpi(string(S.Condition), string(condNames{k}));
end
if ~any(mask), return; end

c = cosd(rotateDeg); s = sind(rotateDeg);
cos0 = S.cosMag(mask);
sin0 = S.sinMag(mask);
S.cosMag(mask) = cos0 .* c - sin0 .* s;
S.sinMag(mask) = sin0 .* c + cos0 .* s;
end
