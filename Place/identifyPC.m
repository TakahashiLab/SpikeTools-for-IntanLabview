%Linearlized version
function [pcList,sinfoO,sinfoS,sRateMap]=identifyPC(ensemble,LTraj,PosT,fstart,varargin)
p = inputParser;
p.addParamValue('percentile', 95, @isnumeric);
p.addParamValue('shufflen', 1000, @isnumeric);
p.addParamValue('event', [], @ismatrix);
p.addParamValue('eventnum', [1 2], @isvector);

p.parse(varargin{:});
Percentile = p.Results.percentile/100;
shuffleN=p.Results.shufflen;
event=p.Results.event;
eventNum=p.Results.eventnum;

fprintf('cell #');

for i=1:size(ensemble,1)
    fprintf('%3d/%3d',i,size(ensemble,1));
    if isempty(event)
        spk=ensemble{i,3};
    else
        spk=extractDelay(ensemble{i,3},event,eventNum);
    end
    
    [oRateMap]=pmap(spk,LTraj,PosT,fstart,'animal','rat');
    
    %%%parpmap
    [sRateMap]=parpmap(spk,LTraj, PosT,fstart,'animal','rat','shuffle',1,'shuffleN',shuffleN);
    %[seq]=pmap(extractDelay(ensemble{i,3},event,nums),LTraj,PosT,0,'animal','rat','verbose',0,'shuffle');

    sinfoO(i)=calcInfo(oRateMap);
    
    for j=1:size(sRateMap,1)
        sinfoS(i,j,:)=calcInfo(squeeze(sRateMap(j,:,:)));
    end
    fprintf('\b\b\b\b\b\b\b');
end
fprintf('\n');
sinfoS=sort(sinfoS,2,'MissingPlacement','first');

pcList=sinfoO>sinfoS(:,size(sinfoS,2)*Percentile)';

return;