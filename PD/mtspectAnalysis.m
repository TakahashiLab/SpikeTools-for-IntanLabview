function [S,f,t]=mtspectAnalysis(lfp,tetnum,Hz,trial)
if nargin==3
  trial=[1 size(lfp,2)];
elseif nargin==2
  trial=[1 size(lfp,2)];
  Hz=25000; 
end

step=floor(Hz/1000);

params.Fs=Hz/step;
params.fpass=[1 40];
params.trialave=0;

data=double(lfp((tetnum-1)*4+1,:));
lfp=decimate(data,step);


%movingwin=[50 5];
movingwin=[50 1];

trial=ceil(trial/step);
loop=length(trial);
duration=min(diff(trial))-1;
%duration=100000;
if loop>2
%  LFP=zeros(duration,loop);
  for i=1:loop
    LFP(:,i)=lfp(trial(i):trial(i)+duration-1);
%    LFP(:,i)=lfp(trial(i)+ceil(duration/2):trial(i)+ceil(duration/2)+10000);
%    LFP(:,i)=lfp(trial(i):trial(i)+10000);
  end
else
  i=1;
  if 1
    LFP=lfp(trial(i):trial(i)+duration-1);
  else
    trialN=10;
    trialD=ceil(size(lfp,2)/trialN)
    for j=1:trialN
      trialD=params.Fs;
      LFP(:,i)=lfp(trial(i)+(j-1)*trialD:trial(i)+j*trialD);
    end
  end
end

%[S,f]=mtspectrumc(LFP,params);
[S,t,f]=mtspecgramc(LFP,movingwin,params);



return;