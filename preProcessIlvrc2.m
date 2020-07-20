%
% preprocessing for dividing ilvrc to electrode, LFP and event data
%
%sleep: 0=all process
%       1=cutting process only
%       2=event process only

function preProcessIlvrc2(basename,sleep,tetrodeNum)
if nargin==1
    sleep=0;
    forceRef=[];
elseif nargin==2
    forceRef=[];
    LEDp=[];
elseif nargin==3
    forceRef=[];
    LEDp=[];
elseif nargin==4
    rf=0;
end


Factor16bit=6553600;
Factor16bit4Event=9930;
uV01=10^7;
uV1=10^6;

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


if sleep==0 | sleep==1 | sleep==2
    
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
            if sleep==0 | sleep==1
                
                [x,ref,sampl]=loadilvrcN(fullfile(dataFolder,d(i).name),tetrodeNum,numOfElectrodes,'tetrode',0,0);
                if ~isempty(forceRef)
                    ref=forceRef;
                end
                           
                e=filterAmp(x,0,sampl);
                e=double(e)./Factor16bit.*uV01;%convert to 0.1 uV
                e=int16(e);
            elseif sleep==2
                if ~exist(savenames{numOfElectrodes+2})
                    event=loadilvrcN(fullfile(dataFolder,d(i).name),tetrodeNum,numOfElectrodes,0,0);
                    event=double(event)./Factor16bit4Event;%convert to V
                end
            end
        end
    end
        
    if sleep==0 | sleep==1
        parsave(savenames{tetrodeNum},'x',e);
        %fprintf('Saving LFPs...\n');
        %parsave(savenames{numOfElectrodes+1},'lfp',lfp,'dlfp',dlfp);
    elseif sleep==2            
        if ~exist(savenames{numOfElectrodes+2})
            fprintf('Saving Events...\n');
            parsave(savenames{numOfElectrodes+2},'event',event);
        end
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