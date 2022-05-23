function [freemem]=unixMem()
    if isunix 
        [r,w] = system('free | grep Mem');
        stats = str2double(regexp(w, '[0-9]*', 'match'));
        freemem = stats(3)/1e6;
    elseif ispc
        [r,w]= memory;
        freemem = w.PhysicalMemory.Total/2^30;
    end
end


