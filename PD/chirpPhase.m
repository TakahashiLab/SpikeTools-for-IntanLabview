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

    for i = 1:length(dInd)
        out = [out -pi:(2 * pi) / dInd(i):pi];
    end
  
    %correction for initial phase
   
    out(1:601)=phase(1:601)*2;
    out(600:900)=-pi:(2*pi)/300:pi;
    
    
end
