%Linearlized version
function [pcList,sinfoO,sinfoS]=identifyPC(ensemble,LTraj,PosT,fstart,varargin)
p = inputParser;
p.addParamValue('percentile', 0.95, @isnumeric);

p.parse(varargin{:});
Percentile = p.Results.percentile;

fprintf('cell #');
for i=1:size(ensemble,1)
    fprintf('%3d/%3d',i,size(ensemble,1));
    [oRateMap]=pmap(ensemble{i,3},LTraj,PosT,fstart,'animal','rat');
    
    %%%parpmap???
    [sRateMap]=parpmap(ensemble{i,3},LTraj, PosT,fstart,'animal','rat','shuffle',1,'shuffleN',100);
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