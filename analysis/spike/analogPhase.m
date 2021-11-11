function xphase=analogPhase(data,th)

[~,locsP]=findpeaks(data,'minpeakheight',th,'minpeakdistance',3,'minpeakwidth',3);
%locs=[1 locs'];

[~,locsN]=findpeaks(-(data-1.5)+1.5,'minpeakheight',th,'minpeakdistance',3,'minpeakwidth',3);

%locsN=[1 locsN'];
xphase=zeros(size(data));

%begining process
if locsP(1)<locsN(1)
  N=locsN(1)-locsP(1);
  c=0;
  for k=locsP(1):-1:1
    xphase(k)=0-pi/N*(c);
    c=c+1;
  end
  
else
  N=locsP(1)-locsN(1);
  xphase(locsN(1):locsP(1))=-pi:pi/N:0;
  c=0;
  for k=locsN(1):-1:1
    xphase(k)=pi-pi/N*(c);
    c=c+1;
  end
  locsN(1)=[];
end


  

if locsN(end)>locsP(end)
  
  %main process
  for i=1:(length(locsP)-1)
    N=locsN(i)-locsP(i);
    xphase(locsP(i):locsN(i))=0:pi/N:pi;
    N=locsP(i+1)-locsN(i);
    xphase(locsN(i):locsP(i+1))=-pi:pi/N:0;
  end
  

%ending process
  N=locsN(end)-locsP(end);
  xphase(locsP(end):locsN(end))=0:pi/N:pi;
  c=0;
  for k=locsN(end):length(xphase)
    xphase(k)=-pi+pi/N*(c);
    c=c+1;
  end
else
  
  %main process
  for i=1:(length(locsP)-1)
    N=locsN(i)-locsP(i);
    xphase(locsP(i):locsN(i))=0:pi/N:pi;
    N=locsP(i+1)-locsN(i);
    xphase(locsN(i):locsP(i+1))=-pi:pi/N:0;
  end
  
  c=0;
  for k=locsP(end):length(xphase)
    xphase(k)=-pi+pi/N*(c);
    c=c+1;
  end
end
  

return;
