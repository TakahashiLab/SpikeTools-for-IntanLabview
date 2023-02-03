function out = chirpPhase(x)
    

    phase = angle(hilbert(x));
    dPhase = diff(phase);

    ind = find(dPhase <- 1);
    ind = [ind; length(x)];

    if ind(1) < 10
        ind(1) = 1;
    else
        ind = [1; ind];
    end

    dInd = diff(ind) - 1;

    out = [-pi];

    for i = 1:(length(dInd))
        out = [out -pi:(2 * pi) / dInd(i):pi];
    end
  
    %correction for start phase
    out(1:15001)=phase(1:15001)*2;
    out(15000:22500)=-pi:(2*pi)/7500:pi;

    out(end-22500:end-15000)=-pi:(2*pi)/7500:pi;
    out(end-15001:end)=phase(end-15001:end)*2;
  
end
