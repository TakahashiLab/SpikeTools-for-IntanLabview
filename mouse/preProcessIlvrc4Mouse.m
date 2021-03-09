%
% preprocessing for dividing ilvrc to electrode, LFP and event data
%
%sleep: 0=all process
%       1=cutting process only 
%       2=event process only
%       3=event process and event cutting
%       4=
%forceRef: When it is not Cell type, number of channel ( Tetrode ...    4-ch1=13)
%forceRef: When it is Cell type, {number of ref channel, number of ...
%    target channels}
%
% ex. split median referencing(left:tetrodes 1-4, right:tetrodes 5-8
%forceRef{1,1}=1;forceRef{1,2}=1:4;forceRef{2,2}=[5:8];
%
% ex. global median referencing(left:tetrodes 1-8)
%forceRef{1,1}=0;
%
%LEDp: a elecrode number which has a optical fiber
%
%ver: 0:mouse old, 1: mouse new2020, 2:
%
function preProcessIlvrc4Mouse(basename,sleep,forceRef,ver,LEDp,rf,gpuFlag)
if nargin==1
  sleep=0;
  forceRef=[];
  rf=0;
  gpuFlag=1;
elseif nargin==2
  forceRef=[];
  LEDp=[];
  rf=0;
  gpuFlag=1;
elseif nargin==3
    ver=0;
    LEDp=[];
    rf=0;
    gpuFlag=1;
elseif nargin==4
    LEDp=[];
    rf=0;
    gpuFlag=1;
elseif nargin==5
    rf=0;
    gpuFlag=1;
elseif nargin==6
    gpuFlag=1;
end


Factor16bit=6553600;
Factor16bit4Event=9930;
uV01=10^7;
uV1=10^6;

numOfMicrowires=4;%tetrode
numOfElectrodes=8;
numOfEvents=5;
Suffix='ilvrc';
SuffixLen=size(Suffix,2)-1;


%load all data
[path,name,ext]=fileparts(basename);

dataFolder=fullfile(path,name);
d=dir(fullfile(path,name));
loop=size(d,1);


SuffixAVI='avi';
SuffixAVILen=size(SuffixAVI,2)-1;
%tracking leds mounted on animal's head
fprintf('tracking leds...\n');

if 1
    %if exist(fullfile(dataFolder,['ledPos.mat']))
  
else
  for i=3:loop
    if strncmp(d(i).name(end-SuffixAVILen:end),SuffixAVI,SuffixAVILen)
      %    Pos=ledtracker(fullfile(dataFolder,d(i).name));
      %    Pos=ledtracker2(fullfile(dataFolder,d(i).name));
      fnp=fullfile(dataFolder,d(i).name);
      Pos=ledtracker3(fnp);
      %load(fullfile(dataFolder,['ledPos.mat']),'Pos');
      %Pos=PosCompensation(Pos,fnp);
      PosOrg=Pos;
      Pos=PosCompensation6(Pos);
      Pos=PosCompT(Pos);
      parsave(fullfile(dataFolder,['ledPos.mat']),'Pos',Pos,'PosOrg',PosOrg);
      clear Pos;
    end 
  end
end

if sleep==0 | sleep==1 | sleep==3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

savenames=cell(1,numOfElectrodes+2);
for i=1:numOfElectrodes
  savenames{i}=fullfile(dataFolder,[name 'e' num2str(i) '.mat']);
end
savenames{numOfElectrodes+1}=fullfile(dataFolder,[name 'LFP.mat']);
xLFP=[];
savenames{numOfElectrodes+2}=fullfile(dataFolder,[name 'Event.mat']);
xEvent=[];
savenames{numOfElectrodes+3}=fullfile(dataFolder,[name 'LFPpower.mat']);


%%loop=4;%for debug
%%while 0
%matlabpool(numOfElectrodes+2);
possibleId=[];
for i=1:loop
    if length(d(i).name)>SuffixLen
        possibleId=[possibleId i];
    end
end
    
fprintf('loading and filtering ilvrc file\n');
for i=possibleId
  if strncmp(d(i).name(end-SuffixLen:end),Suffix,SuffixLen)
    
    if sleep==0 | sleep==1 
      [x,ref,sampl]=loadilvrcN(fullfile(dataFolder,d(i).name),1);
      if ~isempty(forceRef)
	ref=forceRef;
      end

      e=filterAmp(x,0,sampl,gpuFlag);
      x=double(x)./Factor16bit.*uV1;%convert to 1 uV
      e=double(e)./Factor16bit.*uV01;%convert to 0.1 uV
      
      if iscell(ref);
	loop=size(ref,1);
	
	if ref{1,1}==0
	  gmr=median(e,1);%global median referencing 
	  e=e-gmr;
	else% split median referencing
	  for k=1:loop
	    channels=ref{k,2};
	    refNum=ref{k,1};
	    %	    e(channels,:)=e(channels,:)-repmat(e(refNum,:),size(channels,2),1);%referencing 
            pe=e(channels,:);
	    gmr=median(pe,1);
	    e(channels,:)=pe-repmat(gmr,size(channels,2),1);%referencing 
	  end
	end
      else
	e=e-repmat(e(ref,:),size(e,1),1);%referencing 
      end

      e=int16(e);
      
      lfp=filterAmp(double(x),2,sampl,gpuFlag);
      lfp=int16(lfp);
    
      fprintf('convert LFPs\n');
      dlfp=downlfp(lfp,25);%25kHz -> 1kHz
%      dlfp=lfpClean(dlfp);%cleaning using FastICA

      
      if ~exist(savenames{numOfElectrodes+2})
	event=loadilvrcN(fullfile(dataFolder,d(i).name),2);
	event=double(event)./Factor16bit4Event;%convert to V
	%    event=int16(event);
      end
    elseif sleep==3
      if ~exist(savenames{numOfElectrodes+2})
	event=loadilvrcN(fullfile(dataFolder,d(i).name),2);
	event=double(event)./Factor16bit4Event;%convert to V
	%    event=int16(event);
      end
    end
  end
end

if sleep==0 | sleep==1 | sleep==3
  %parfor j=1:(numOfElectrodes+2)
  
  if sleep==0 | sleep==1
    for j=1:(numOfElectrodes)
      fprintf('Saving electrode # %d...\n',j);
      x=e(numOfMicrowires*(j-1)+1:numOfMicrowires*j,:);
      parsave(savenames{j},'x',x); 
    end	

    fprintf('Saving LFPs...\n');
    %x=lfp;
    parsave(savenames{numOfElectrodes+1},'lfp',lfp,'dlfp',dlfp); 
  end
  
  %x=event;
  if ~exist(savenames{numOfElectrodes+2})
    fprintf('Saving Events...\n');
    parsave(savenames{numOfElectrodes+2},'event',event); 
  end
end
%matlabpool close;
%return;
end


%if sleep==0 | sleep==2 | sleep==3
if 1
  fprintf('Saving Position and Task...\n');
  %TaskEvents=extractTask(basename,fullfile(dataFolder,[name 'Event.mat']));
  %save(fullfile(dataFolder,[name 'task.mat']),'TaskEvents'); 
  %  if exist(fullfile(dataFolder,['ledPos.mat']))
  %makePositionIntan(basename);
  if exist(fullfile(dataFolder,['DLC.mat']))
      load DLC.mat Traj;
      %PosT=extractPos(fn,Traj,4);
      [PosT,Traj]=extractPos(basename,Traj,3,ver);
      sfn=fullfile(dataFolder,'positions.mat');
      Pos=Traj;
      save(sfn,'Pos','PosT');
  end
  
  if rf==1
    [~,TrialT]=extractTaskIntanlv(basename,rf);%produce timing.mat
                                               %load(fullfile(dataFolder,['positions.mat']),'Pos','PosT');
    segPara=extractRFevent(TrialT,Pos,PosT);
    sfn=fullfile(dataFolder,'timing.mat');
    save(sfn,'segPara','TrialT');
  else
    extractTaskIntanlv(basename);%produce timing.mat
  end
end

if sleep~=3
  %LFP power
  load(fullfile(dataFolder,['positions.mat']));
  if isempty(LEDp)
    mP=LFPpower(dlfp,PosT);
  else
    ass{1}=LEDp;
    if LEDp<5
      ass{2}=setdiff(1:4,LEDp);
      ass{3}=5:8;
    else
      ass{2}=setdiff(5:8,LEDp);
      ass{3}=1:4;
    end
    mP=LFPpower(dlfp,PosT,ass);
  end
  parsave(savenames{numOfElectrodes+3},'mP',mP);
end



return;

