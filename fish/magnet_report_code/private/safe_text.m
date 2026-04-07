function txt = safe_text(x)
if isstring(x), x = x(1); end
if iscell(x),   x = x{1}; end
if isnumeric(x), x = num2str(x); end
txt = char(x);
txt = strrep(txt, char(8211), '-');
txt = strrep(txt, char(8212), '-');
end
