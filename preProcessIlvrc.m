%
% preprocessing for dividing ilvrc to electrode, LFP and event data
%
%sleep: 0=all process
%       1=cutting process only
%       2=event process only
%       3=event process and event cutting
%       4=kilosort process only
%forceRef: When it is not Cell type, number of channel ( Tetrode ...    4-ch1=13)
%forceRef: When it is Cell type, {number of ref channel, number of ...
%    target channels}
% ex.) forceRef{1,1}=1;
%      forceRef{1,2}=5:8;
%LEDp: a elecrode number which has a optical fiber
function preProcessIlvrc(basename,sleep,forceRef,LEDp,rf)
if nargin==1
    sleep=0;
    forceRef=[];
elseif nargin==2
    forceRef=[];
    LEDp=[];
elseif nargin==3
    LEDp=[];
elseif nargin==4
    rf=0;
end


Factor16bit=6553600;
Factor16bit4Event=9930;
uV01=10^7;

numOfMicrowires=4;%tetrode
numOfElectrodes=16;%64 channels
numOfEvents=5;
Suffix='ilvrc';
SuffixLen=size(Suffix,2)-1;


%load all data
[path,name,ext]=fileparts(basename);
dataFolder=fullfile(path,name);
d=dir(fullfile(path,name));
loop=size(d,1);


if sleep==0 | sleep==1 | sleep==3 | sleep==4%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    savenames=cell(1,numOfElectrodes+2);
    for i=1:numOfElectrodes
        savenames{i}=fullfile(dataFolder,[name 'e' num2str(i) '.mat']);
    end
    savenames{numOfElectrodes+1}=fullfile(dataFolder,[name 'LFP.mat']);
    xLFP=[];
    savenames{numOfElectrodes+2}=fullfile(dataFolder,[name 'Event.mat']);
    xEvent=[];
    savenames{numOfElectrodes+3}=fullfile(dataFolder,[name 'LFPpower.mat']);
    
    savenames{numOfElectrodes+4}=fullfile(dataFolder,[name '.bin']);
    
    fprintf('loading and filtering ilvrc file\n');
    possibleId=[];
    for i=1:loop
        if length(d(i).name)>SuffixLen
            possibleId=[possibleId i];
        end
    end
    
    
    for i=possibleId
        if strncmp(d(i).name(end-SuffixLen:end),Suffix,SuffixLen)
            if sleep==4
                [x,ref,sampl]=loadilvrcN(fullfile(dataFolder,d(i).name),1,numOfElectrodes,'tetrode',0,1);
                if ~exist(savenames{numOfElectrodes+2})
                    event=loadilvrcN(fullfile(dataFolder,d(i).name),2);
                    event=double(event)./Factor16bit4Event;%convert to V
                end
            elseif sleep==0 | sleep==1
                
                [x,ref,sampl]=loadilvrcN(fullfile(dataFolder,d(i).name),1,numOfElectrodes,'tetrode',0,1);
                if ~isempty(forceRef)
                    ref=forceRef;
                end
                
                
                e=filterAmp(x,0,sampl);
                e=double(e)./Factor16bit.*uV01;%convert to 0.1 uV
                
                if iscell(ref);
                    
                    if ref{1,1}==0
                        gmr=median(e,1);%global median referencing
                        e=e-gmr;
                    else% split median referencing
                        loop1=size(ref,1);
                        for k=1:loop1
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
                
                lfp=filterAmp(double(x),2,sampl);
                lfp=int16(lfp);
                
                fprintf('convert LFPs\n');
                dlfp=downlfp(lfp,25);%25kHz -> 1kHz
                %      dlfp=lfpClean(dlfp);%cleaning using FastICA
                
                if ~exist(savenames{numOfElectrodes+2})
                    event=loadilvrcN(fullfile(dataFolder,d(i).name),2);
                    event=double(event)./Factor16bit4Event;%convert to V
                end
            elseif sleep==3
                if ~exist(savenames{numOfElectrodes+2})
                    event=loadilvrcN(fullfile(dataFolder,d(i).name),2);
                    event=double(event)./Factor16bit4Event;%convert to V
                end
            end
        end
    end
    
    if sleep==0 | sleep==1 | sleep==3 | sleep==4
        
        
        if sleep==0 | sleep==1
            saveKilosort(savenames{numOfElectrodes+4},x);
            for j=1:(numOfElectrodes)
                fprintf('Saving electrode # %d...\n',j);
                x=e(numOfMicrowires*(j-1)+1:numOfMicrowires*j,:);
                parsave(savenames{j},'x',x);
            end
            fprintf('Saving LFPs...\n');
            parsave(savenames{numOfElectrodes+1},'lfp',lfp,'dlfp',dlfp);
        elseif sleep==4
            saveKilosort(savenames{numOfElectrodes+4},x);
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

fprintf('Saving Position and Task...\n');

if 0
    if rf==1
        [~,TrialT]=extractTaskIntanlv(basename,rf);%produce timing.mat
        load(fullfile(dataFolder,['positions.mat']),'Pos','PosT');
        segPara=extractRFevent(TrialT,Pos,PosT);
        sfn=fullfile(dataFolder,'timing.mat');
        save(sfn,'segPara','TrialT');
    else
        extractTaskIntanlv(basename);%produce timing.mat
    end
end


return;
%%%%%%%%%
function saveKilosort(fn,x)
fid=fopen(fn,'w');
x=int16(double(x)/10);
fwrite(fid,x,'int16');
fclose(fid);
return;