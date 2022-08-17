%dn: 2018,2019,2020
function [phaseHistPyr,phaseHistInt,PyrIntList,PyrIntListStim,FRs,TPs,SWs,PIs]=contPD(dn,an,varargin)

p = inputParser;
p.addParamValue('method', 'gainmap', @ischar);
p.addParamValue('wavetype', 'chirp', @ischar);
p.addParamValue('intensity', 5, @isnumeric);
p.addParamValue('cellclass', 0, @isnumeric);%cell classification mode


p.parse(varargin{:});
method = p.Results.method;
waveType = p.Results.wavetype;
Intensity=p.Results.intensity;
cellClass=p.Results.cellclass;


alpha=0.05;
[homePath,dataPath]=PDdataC(dn,an);
TimingData='timing.mat';
EnsembleData='ensemble.mat';
details=dataPath;

loop=size(dataPath,1);
params=cell(loop,2);

phaseHistPyr=[];
phaseHistInt=[];
offSetPyr=0;
offSetInt=0;

normalPyr=[];
normalInt=[];
PDpyr=[];
PDint=[];

normalPyrStim=[];
normalIntStim=[];
PDpyrStim=[];
PDintStim=[];

PyrIntList=[];
PyrIntListStim=[];

FRs=[];
TPs=[];
SWs=[];
PIs=[];




for i=1:loop

    if dataPath{i,7}==Intensity
        fprintf('analyzing %s\n',dataPath{i,1});
        [~,name]=fileparts(dataPath{i,1});

        switch lower(method)
          case 'behavior',
            loadname=fullfile(homePath,dataPath{i,1},TimingData);
            load(loadname,'segPara','TrialT');
            if dataPath{i,4}<=4%LFP reference was on the left hemisphere
                if dataPath{i,4}==dataPath{i,3}
                    LEDout=[LEDout;cellfun(@median,segPara(:,2))'];%1:speed, 2:maxdist     
                else
                    LFPoutL=[LFPoutL;cellfun(@median,segPara(:,2))'];%1:speed, 2:maxdist
                end
            else%Optical fiber was implanted in the right hemisphere
                if dataPath{i,4}==dataPath{i,3}
                    LEDout=[LEDout;cellfun(@median,segPara(:,2))'];%1:speed, 2:maxdist  
                else
                    LFPoutR=[LFPoutR;cellfun(@median,segPara(:,2))'];%1:speed, 2:maxdist     
                end
            end
          case 'gainmap',

            loadname=fullfile(homePath,dataPath{i,1},TimingData);
            load(loadname,'preSilentTime','postSilentTime','chirpSeqTime','noiseSeqTime');
            loadname=fullfile(homePath,dataPath{i,1},EnsembleData);
            load(loadname,'ensemble','an','en');
            %classify pyr and int
            [pyr,interneuron,fr,tp,sw]=checkFR(ensemble);
            fprintf('real interneuron\n');
            interneuron=union(dataPath{i,8},dataPath{i,10});%%%%%CC+tag
            pyrTag=setdiff(dataPath{i,11},interneuron);%%conflict
            pyr=union(dataPath{i,9},pyrTag);%%%CC+tag
                                  
            switch  lower(waveType)
              case 'chirp',
                SeqTime=chirpSeqTime;
              case 'noise',
                SeqTime=noiseSeqTime;
            end
            
            if cellClass
                FRs=[FRs; fr];
                TPs=[TPs; tp];
                SWs=[SWs; sw];
                PIs=[PIs; sparse(dataPath{i,8},1,1,length(ensemble),1)];
            else
                loadname=fullfile(homePath,dataPath{i,1},[name 'Event.mat']);
                load(loadname,'event');


                %calc gain map
                [pPyr,pInt]=batchGainMap(event,ensemble(:,3),SeqTime,1,[preSilentTime;postSilentTime],pyr,interneuron);

                %cell classification
                normal=find(tetrodeMap(an,en)<=4);
                if ~isempty(pyr)
                    %pyramidal
                    
                    [~,ind]=intersect(pyr,1:normal(end));
                    normalPyr=[normalPyr; ind+offSetPyr];
                    if dataPath{i,3}<=4
                        normalPyrStim=[normalPyrStim; ind+offSetPyr];
                    end

                    [~,ind]=setdiff(pyr,1:normal(end));
                    PDpyr=[PDpyr; ind+offSetPyr];
                    if dataPath{i,3}>=5
                        PDpyrStim=[PDpyrStim; ind+offSetPyr];
                    end
                end

                %interneuron
                [~,ind]=intersect(interneuron,1:normal(end));
                normalInt=[normalInt; ind+offSetInt];
                if dataPath{i,3}<=4
                    normalIntStim=[normalIntStim; ind+offSetInt];
                end

                [~,ind]=setdiff(interneuron,1:normal(end));
                PDint=[PDint; ind+offSetInt];
                if dataPath{i,3}>=5
                    PDintStim=[PDintStim; ind+offSetInt];
                end

                offSetPyr=offSetPyr+length(pyr);
                offSetInt=offSetInt+length(interneuron);
                
                phaseHistPyr=cat(3,phaseHistPyr,pPyr);
                phaseHistInt=cat(3,phaseHistInt,pInt);
            end

        end    
    end

end


if ~cellClass
    PyrIntList{1}=normalPyr;
    PyrIntList{2}=normalInt;
    PyrIntList{3}=PDpyr;
    PyrIntList{4}=PDint;
    PyrIntListStim{1}=normalPyrStim;
    PyrIntListStim{2}=normalIntStim;
    PyrIntListStim{3}=PDpyrStim;
    PyrIntListStim{4}=PDintStim;
end

return;