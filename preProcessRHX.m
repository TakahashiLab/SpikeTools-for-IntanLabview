%
% preprocessing for dividing ilvrc to electrode, LFP and event data
%
%sleep: 0=all process
%       1=cutting process only
%       2=event process only
%       3=event process and event cutting
%forceRef: When it is not Cell type, number of channel ( Tetrode ...    4-ch1=13)
%forceRef: When it is Cell type, {number of ref channel, number of ...
%    target channels}
% ex.) forceRef{1,1}=1;
%      forceRef{1,2}=5:8;
%LEDp: a elecrode number which has a optical fiber
function preProcessRHX(basename,sleep,forceRef,gpuFlag)
if nargin==1
    sleep=0;
    forceRef=[];
    gpuFlag=1;
elseif nargin==2
    forceRef=[];
    gpuFlag=1;
elseif nargin==3
    gpuFlag=1;
end


Factor16bit=6553600;
Factor16bit4Event=9930;
uV01=10^7;
uV1=10^6;

sampl=25000;%25kHz sampling

numOfMicrowires=4;%tetrode
numOfElectrodes=16;%64 channels
numOfEvents=5;
Suffix='rhd';
SuffixLen=size(Suffix,2)-1;


%load all data
[path,name,ext]=fileparts(basename);
dataFolder=fullfile(path,name);
d=dir(fullfile(path,name));
loop=size(d,1);

possibleId=[];
for i=1:loop
    if length(d(i).name)>SuffixLen
        possibleId=[possibleId i];
    end
end
    
for i=possibleId
    if strncmp(d(i).name(end-SuffixLen:end),Suffix,SuffixLen)
        filename=fullfile(dataFolder,d(i).name);
    end
end
            

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
    
    fprintf('loading RHX file:%s\n',filename);
    
    if sleep==0 | sleep==1
        [x,event,t]=loadRHD(filename);
        if ~isempty(forceRef)
            ref=forceRef;
        end
        e=filterAmp(x,0,sampl,gpuFlag);          
        %%x=double(x)./Factor16bit.*uV1;%convert to 1 uV

        lfp=filterAmp(double(x),2,sampl,gpuFlag);
        lfp=int16(lfp);
        fprintf('convert LFPs\n');
        dlfp=downlfp(lfp,25);%25kHz -> 1kHz
        fprintf('Saving LFPs...\n');
        parsave(savenames{numOfElectrodes+1},'lfp',lfp, ...
                'dlfp',dlfp);
        clear x lfp dlfp;
                
        %e=double(e)./Factor16bit.*uV01;%convert to 0.1 uV
                
        if iscell(ref);
            loop=size(ref,1);
            if ref{1,1}==0
                gmr=median(e,1);%global median referencing
                e=e-gmr;
            else% split median referencing
                for k=1:loop
                    channels=ref{k,2};
                    refNum=ref{k,1};
                    pe=e(channels,:);
                    gmr=median(pe,1);
                    e(channels,:)=pe-repmat(gmr,size(channels,2),1); ...
                    %referencing
                end
            end
        else
            e=e-repmat(e(ref,:),size(e,1),1);%referencing
        end
                
        e=int16(e);
                
                
        if ~exist(savenames{numOfElectrodes+2})
            %%event=loadilvrcN(fullfile(dataFolder,d(i).name),2);
            %%event=double(event)./Factor16bit4Event;%convert to V
        end
    elseif sleep==3
        if ~exist(savenames{numOfElectrodes+2})
            %% event=loadilvrcN(fullfile(dataFolder,d(i).name),2);
            %%event=double(event)./Factor16bit4Event;%convert to V
        end
    end

    
    if sleep==0 | sleep==1 | sleep==3 | sleep==4
        
        
        if sleep==0 | sleep==1

            for j=1:(numOfElectrodes)
                fprintf('Saving electrode # %d...\n',j);
                x=e(numOfMicrowires*(j-1)+1:numOfMicrowires*j,:);
                parsave(savenames{j},'x',x);
            end

        end
        
        %x=event;
        if ~exist(savenames{numOfElectrodes+2})
            fprintf('Saving Events...\n');
            parsave(savenames{numOfElectrodes+2},'event',event,'t',t);
        end
    end
    %matlabpool close;
    %return;
end


if exist(fullfile(dataFolder,['DLC.mat']))
    if sleep==3
        load(savenames{numOfElectrodes+2},'event','t');
    end
    load DLC.mat Pos;
    Pos=Pos';
    [NosePork,Treadmill,Pos,PosT]=extractRHX(event,t,Pos);   
    PosT=PosT';
    save event.mat NosePork Treadmill Pos PosT;
    sfn=fullfile(dataFolder,'positions.mat');
    save(sfn,'Pos','PosT');
end

if 0
fprintf('Saving Position and Task...\n');

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
x=int16(x);
fwrite(fid,x,'int16');
fclose(fid);
return;