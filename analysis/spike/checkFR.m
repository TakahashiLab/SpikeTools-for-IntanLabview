%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pyr,interneuron,FRs,TPs,SWs]=checkFR(Out)
Th=5;
loop=size(Out,1);

if size(Out,2)>1
  inputType=1;
else
  inputType=0;
end
FRs=[];
TPs=[];
SWs=[];
interneuron=[];
kHz=25;

E=realmin;
S=realmax;
for i=1:loop
    if inputType
    tmps=Out{i,3};
  else
    tmps=Out{i};
  end
  if ~isempty(tmps)
    if tmps(end) > E
      E=tmps(end);
    end
    if tmps(1)<S
      S=tmps(1);
    end
  end
end  

len=(E-S)/(kHz*1000);

for i=1:loop
  if inputType
    tmps=Out{i,3};
    wf=Out{i,1};
  else
    tmps=Out{i};
  end
    
  if ~isempty(tmps)

    FR=size(tmps,2)/double(len);

    %if FR >Th
    %fprintf('deleting #%d cell %2.1fHz\n',i,FR);
    %interneuron=[interneuron i];
    %    end
    %trough-peak




    len=length(tmps);



    wf=reshape(wf,4,size(wf,2)/len,len);
    wf=mean(wf,3);
    [~,wid]=max(abs(max(abs(wf),[],2)));

    thiswave=wf(wid,:);
    [minval,minpos] = min(thiswave);
    minpos = minpos(1);
    [maxval,maxpos] = max(thiswave);
    [dummy,maxpos] = max(thiswave(minpos+1:end));
    if isempty(maxpos)
        warning('Your Waveform may be erroneous')
        maxpos = 1
    end
    maxpos=maxpos(1);
    maxpos = maxpos+minpos;
    TP = maxpos-minpos; %In number of samples


    w = wf(wid,:)';
    w = [w(1)*ones(1000,1);w;w(end)*ones(1000,1)];
    [wave f t] = getWavelet(w,20000,500,3000,128);
    %We consider only the central portion of the wavelet because we
    %haven't filtered it before hand (e.g. with a Hanning window)
    wave = wave(:,int16(length(t)/4):3*int16(length(t)/4));
    %Where is the max frequency?
    [maxPow ix] = max(wave);
    [dumy mix] = max(maxPow);
    ix = ix(mix);
    SW = 1000/f(ix);

    
  else
    FR=[];
    TP=[];
    SW=[];
  end
  
  FRs=[FRs; FR];
  TPs=[TPs; TP];
  SWs=[SWs; SW];
end

TPs=TPs/kHz;

x=TPs;
y=SWs;

xx = [0 0.8];
yy = [0.8 0];
m = diff( yy ) / diff( xx );
b = yy( 1 ) - m * xx( 1 );  % y = ax+b
pyr = y>= m*x+b;
interneuron = ~pyr;

pyr=find(pyr)';
interneuron=find(interneuron)';

%interneuron=find(TPs<0.4)';
%pyr=setdiff(1:loop,interneuron);

return;


