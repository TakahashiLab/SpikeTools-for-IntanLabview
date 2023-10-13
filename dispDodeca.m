%
%function dispDodeca(Out,num,varargin)
%dispDodeca(Out,1,'mode',1)
function [MF,width,FR]=dispDodeca(Out,num,varargin)

Ts=[];

kHz=25;
step=size(Out{num,1},2)/size(Out{num,3},2);
%step=40;
limitLen=1000;
%AmpFactor=2500/(2^15);old
AmpFactor=1;
verbose=1;
average=0;

modeP=2;%dodecatrode mode
raw=0;%Independent component
silent=1;

for i=1:2:(length(varargin)-1)
  % change the value of parameter
  switch lower(varargin{i})
    case 'khz'
      kHz=varargin{i+1};
    case 'limitlen'
      limitLen=varargin{i+1};
    case 'mode'
      raw=varargin{i+1};
    case 'timerange'
      Ts=varargin{i+1};
      Ts=Ts*60*1000*20;
    case 'ampfactor'
      %      AmpFactor=varargin{i+1}/(2^15);
      AmpFactor=varargin{i+1};
    case 'silent'
      silent=varargin{i+1};
    case 'average'
      average=varargin{i+1};
    case 'verbose'
      verbose=varargin{i+1};
    otherwise
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized parameter: ''' varargin{i} '''']);
      return;
  end
end


if raw
  tmp=Out{num,4};  
else
  tmp=Out{num,1};
end

tmps=Out{num,3};


MF=makeMeanWF(double(tmp),step);
[dum,targetId]=min(min(MF,[],2));
[dum,maxId]=max(MF(targetId,:));
[dum,minId]=min(MF(targetId,:));

if silent
  if ~isempty(Ts)
    seq=find(tmps>Ts(1) & tmps<Ts(2));
    wf=tmp(:,1+(seq(1)-1)*step:seq(end)*step);
    core(wf,limitLen,step,kHz,AmpFactor);
  else
      if verbose
          if average
              core(tmp,limitLen,step,kHz,AmpFactor,MF);
          else
              core(tmp,limitLen,step,kHz,AmpFactor);
          end
      end
  end

  if ~isempty(Ts)
    seq=find(tmps>Ts(2) & tmps<Ts(3));
    wf=tmp(:,1+(seq(1)-1)*step:seq(end)*step);
    figure;
    core(wf,limitLen,step,kHz,AmpFactor);
  end
end



if silent

width=abs(maxId-minId)/kHz;
fprintf('spikeWidth=%1.2fmsec\n',width);
len=((tmps(end)-tmps(1))/kHz/1000);

FR=size(tmps,2)/double(len);
fprintf('Firing Rate=%1.2fHz\n',FR);

if 0
    %if exist('ranksum') & exist('ttest2')
  %[out,num]=getBurst(tmp,tmps,step,kHz);
  [p1,dum,dum,p]=getBurst(tmp,tmps,step,kHz);

  if p>0 & p <= 0.05 
    fprintf('Dendrite:%f\n',p);
  elseif p < 1
    fprintf('Soma:%f\n',p);
  else
    fprintf('Unidentifed\n');
  end
end
end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function core(tmp,limitLen,step,kHz,AmpFactor,MF)
if nargin==5
    MF=[];
end
    
[maxstep,Y]=size(tmp);

startLoop=1;
if length(limitLen)==2
  startLoop=limitLen(1);
  endLoop=Y/step;
  if endLoop> limitLen
    endLoop=limitLen(2);
  end
else
  endLoop=Y/step;
  if endLoop> limitLen
    endLoop=limitLen;
  end
end


MaxAmp=max(tmp(:)*AmpFactor*1.01);
MinAmp=min(tmp(:)*AmpFactor*1.01);

clf;
for loop=startLoop:endLoop
  hold on ;
  for dit = 1:maxstep
    xrange = ((step+2) * (dit-1)) + (1:step);
    range=((loop-1)*step+1) : (step*loop);

    plot(xrange, tmp(dit,range)*AmpFactor);
  end
end

if ~isempty(MF)
    for dit = 1:maxstep
        xrange = ((step+2) * (dit-1)) + (1:step);
        range=((step) * (dit-1)) + 1:((step) * (dit-1)) + step;
        plot(xrange, MF(range)*AmpFactor,'LineWidth',10);
    end      
end


axis([0 (step+2)*(maxstep)-2 MinAmp MaxAmp]);
axis off;

hold on;
plot([(step+2)*(maxstep-1)-2 (step+2)*(maxstep-1)-2+kHz],[MinAmp MinAmp],'k');
plot([(step+2)*(maxstep-1)-2 (step+2)*(maxstep-1)-2],[MinAmp MinAmp+10],'k');

X=(step+2)*(maxstep-1)-2;
Y=double(MinAmp)*1.1;
text(X,Y,'1ms');

Y=double(MinAmp)*0.9;
text(X,Y,'10uV');
return;
