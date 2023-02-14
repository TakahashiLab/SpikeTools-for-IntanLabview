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
  
    %interval 60 sec 1kHz
    %correction for start phase
    interval=60*1000;
    loop=floor(length(out)/interval);

    for q=1:loop
      tmp=out(1+(q-1)*interval:q*interval);
      
      tmp(1:601)=[0:pi/600:pi];
      tmp(interval-600:interval)=[-pi:pi/600:0];

      out(1+(q-1)*interval:q*interval)=tmp;
    end

end
