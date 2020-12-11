function [PulseT,TrialT]=extractTaskIntanlv(filename,rf)

if nargin==1
  rf=0;
else
  PulseT=[];
end

Hz=25000;%
jitter=300000;
%Silent duration 600sec
silentDuration=Hz*600;

%pulse 7 patterns, up and down, 10sec duration(up/down) 
pulseDuration=Hz*7*(10*2);

%chirp pattern, up and down, 30 sec duraion(up/down)
chirpDuration=Hz*(30*2);

%noise pattern, 10 sec duration
noiseDuration=Hz*10;

%tagging pattern, 60sec duration
tagDuration=Hz*60;

[path,name,ext]=fileparts(filename);
Efilename = fullfile(path,[name '/' name 'Event.mat']);

if ~isempty(who('-file',Efilename,'event'))
  load(Efilename,'event');
else
  load(Efilename,'x');
  event=x;
end
%Th=2^14;
Pulse=event(1,:);
Trial=event(2,:);
Th=max(Pulse(:))/2;


[TrialT]=getTimes(Trial,Th);


if rf~=1
  
dTr=diff(TrialT);

%pulse
pulseSeq=find(dTr>pulseDuration-jitter & dTr<pulseDuration+jitter);
pulseStartTime=TrialT(pulseSeq(1));
pulseStopTime=TrialT(pulseSeq(end)+1);
pulseSeqTime=TrialT(pulseSeq);
fprintf('# of trials for 7 pulse patterns=%d\n',length(pulseSeq));
pulseTime=[pulseStartTime pulseStopTime];

%chirp
chirpSeq=find(dTr>chirpDuration-jitter & dTr<chirpDuration+jitter);
chirpSeq(find(diff(chirpSeq)>1)+1)=[];
chirpSeqTime=TrialT(chirpSeq);
fprintf('# of trials for chirp patterns=%d\n',length(chirpSeq));
chirpStartTime=TrialT(chirpSeq(1));
chirpStopTime=TrialT(chirpSeq(end)+1);
chirpTime=[chirpStartTime chirpStopTime];

%noise
noiseSeq=find(dTr>noiseDuration-jitter & dTr<noiseDuration+jitter);
noiseSeq(find(diff(noiseSeq)>1)+1)=[];
noiseSeqTime=TrialT([noiseSeq noiseSeq(end)+1]);%modified
fprintf('# of trials for noise patterns=%d\n',length(noiseSeqTime));
noiseStartTime=TrialT(noiseSeq(1));
noiseStopTime=TrialT(noiseSeq(end)+1)+noiseDuration;
noiseTime=[noiseStartTime noiseStopTime];


%tagging
tagSeq=find(dTr>tagDuration-jitter & dTr<tagDuration+jitter);
tagSeq=tagSeq(end);
tagStartTime=TrialT(tagSeq);
tagStopTime=TrialT(tagSeq+1);
tagTime=[tagStartTime tagStopTime];

%pre Silent 
preSilentStopTime=TrialT(pulseSeq(1));
preSilentStartTime=TrialT(pulseSeq(1))-silentDuration;
preSilentTime=[preSilentStartTime preSilentStopTime];

%post Silent 
%postSilentStopTime=TrialT(noiseSeq(end)+1)+silentDuration;
%postSilentStartTime=TrialT(noiseSeq(end)+1);
postSilentStopTime=noiseStopTime+silentDuration;
postSilentStartTime=noiseStopTime+1;
postSilentTime=[postSilentStartTime postSilentStopTime];

[PulseT]=getTimes(Pulse,Th);

pulsePoint=PulseT(find(PulseT>pulseTime(1) & PulseT<pulseTime(2)));
tagPoint=PulseT(find(PulseT>tagTime(1) & PulseT<tagTime(2)));

end


dataFolder=fullfile(path,name);
sfn=fullfile(dataFolder,'timing.mat');

if rf~=1
  save(sfn,'preSilentTime','postSilentTime','noiseTime','chirpTime','pulseTime','tagPoint','pulsePoint','tagTime','chirpSeqTime','pulseSeqTime','noiseSeqTime');

end
return;
%LED=x(3,:);
%Th=max(LED(:))/2;

%[PosT]=getTimes2(LED,Th);
%sfn=fullfile(dataFolder,'ledPos.mat');
%save(sfn,'PosT','-append');

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%get times
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PosT=getTimes2(x,Th)

P=find(x>Th);
P=[0 P];

tmp=find((diff(P)~=1));
tmp=unique([tmp tmp+1]);
tmp(1)=[];

PosT=P(tmp);

%remove oversampled data(interval <5msec)
ro=find(diff(PosT)<125);%125
ro=unique([ro ro+1]);

if size(ro,2)>1
  fprintf('removing oversampled data:%d...\n',size(ro,2));
end
PosT(ro)=[];

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%get times
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PosT=getTimes(x,Th)

P=find(x>Th);
P=[0 P];

tmp=find((diff(P)~=1));
tmp=unique([tmp tmp+1]);
tmp(1)=[];

PosT=P(tmp);

%remove oversampled data(interval <5msec)
ro=find(diff(PosT)<=1255);%125
ro=unique([ro+1]);

%if size(ro,2)>1
%  fprintf('removing oversampled data:%d...\n',size(ro,2));
%end
PosT(ro)=[];

return;
