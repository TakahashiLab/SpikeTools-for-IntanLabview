function out=multiCFC(dlfp,EpsEdge,elecN)
if nargin<3
  elecN=8;
end

numOfElec=4;
if ~isempty(EpsEdge)
  EpsLen=size(EpsEdge,2);
else
  EpsLen=1;
end

for Eps=1:size(EpsEdge,2)
  figure;
  for i=1:elecN
    subplot(4,2,i)
    
    if ~isempty(EpsEdge)
      cEps=floor(EpsEdge{Eps}./25);%%%
      duration=[];
      for j=1:size(cEps,1)
	duration=[duration cEps(j,1):cEps(j,2)];
      end
    else
      duration=1:size(dlfp,2);
    end
    [xphase,xamp]=PAC(dlfp(1+(i-1)*numOfElec,duration));
    out=getCFC(xphase(1:20,:)',xamp');
    imagesc(out.MI);
  end
end

return;