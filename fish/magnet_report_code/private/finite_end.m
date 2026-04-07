function x = finite_end(ws, wl)
if isfinite(wl), x = ws+wl; else, x = Inf; end
end
