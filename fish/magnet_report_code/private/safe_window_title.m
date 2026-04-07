function txt = safe_window_title(prefix, ws, wl)
txt = safe_text(sprintf('%s (Window %gs-%gs)', prefix, ws, finite_end(ws, wl)));
end
