function [seq,order,oc_map,LTraj,LTrajS,TMpoints]=SequencePmap(ensemble,event,nums,Traj,PosT,ThS,varargin)

p = inputParser;
p.addParamValue('outputmode', 'raw', @ischar);
p.addParamValue('ltraj', [], @ismatrix);
p.addParamValue('jitterpre', 0, @isnumeric);
p.addParamValue('jitterpost', 0, @isnumeric);
p.addParamValue('shuffle', 0, @isnumeric);
p.parse(varargin{:});
outputMode = p.Results.outputmode;
LTraj = p.Results.ltraj;
jitterPre=p.Results.jitterpre;
jitterPost=p.Results.jitterpost;
shuffle=p.Results.shuffle;

if nargin<=5
    ThS=2.5;
end

if isempty(LTraj)
    [LTraj,LTrajS,TMpoints]=fitLinear(Traj);
else
    LTrajS=[];
    TMpoints=[];
end

loop=size(ensemble,1);
for i=1:loop
    fprintf('%d/%d\n',i,loop);
    if isempty(nums)
        if shuffle
            [seq(i,:,:)]=pmap(ensemble{i,3},LTraj,PosT,0,'animal','rat','speed',ThS,'verbose',0,'shuffle',shuffle);
            oc_map=[];
        else
            [seq(i,:),~,~,oc_map]=pmap(ensemble{i,3},LTraj,PosT,0,'animal','rat','speed',ThS,'verbose',0);
        end
    else
        if shuffle
            [seq(i,:,:)]=pmap(ensemble{i,3},LTraj,PosT,0,'animal','rat','speed',ThS,'verbose',0,'shuffle',shuffle);
            oc_map=[];
        else
            [seq(i,:),~,~,oc_map]=pmap(extractDelay(ensemble{i,3},event,nums,'jitterPre',jitterPre,'jitterPost',jitterPost),LTraj,PosT,0,'animal','rat', ...
                 'speed',ThS,'verbose',0);
        end
    end
end

[maxSeq,ind]=max(seq,[],2);
[~,ind]=sort(ind);

switch (outputMode)
  case 'normalized',
    seq=seq./maxSeq;
  case 'raw',
end

order=ind;

return;