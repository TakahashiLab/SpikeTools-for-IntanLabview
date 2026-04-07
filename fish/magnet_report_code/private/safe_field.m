function s = safe_field(x)
s = char(safe_text(x)); s(~isstrprop(s,'alphanum')) = '_';
end
