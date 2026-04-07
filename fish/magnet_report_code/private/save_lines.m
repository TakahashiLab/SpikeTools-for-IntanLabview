function save_lines(path, lines_str)
try
    writelines(lines_str, path);
catch
    fid = fopen(path, 'w');
    if fid<0, warning('Cannot open %s', path); return; end
    for i=1:numel(lines_str)
        fprintf(fid, '%s\n', safe_text(lines_str(i)));
    end
    fclose(fid);
end
end
