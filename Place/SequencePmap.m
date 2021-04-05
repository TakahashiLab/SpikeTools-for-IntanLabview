function [seq,order,oc_map,LTraj,LTrajS,TMpoints]=SequencePmap(ensemble,event,nums,Traj,PosT,ThS,varargin)

p = inputParser;
p.addParamValue('outputmode', 'raw', @ischar);
p.addParamValue('ltraj', [], @ismatrix);
p.addParamValue('jitterpre', 0, @isnumeric);
p.addParamValue('jitterpost', 0, @isnumeric);
p.parse(varargin{:});
outputMode = p.Results.outputmode;
LTraj = p.Results.ltraj;
jitterPre=p.Results.jitterpre;
jitterPost=p.Results.jitterpost;

if nargin<=5
    ThS=2.5;
    if isempty(LTraj)
        [LTraj,LTrajS,TMpoints]=fitLinear(Traj);
    end
elseif nargin<=6
    if isempty(LTraj)
        [LTraj,LTrajS,TMpoints]=fitLinear(Traj);
    end
end


loop=size(ensemble,1);
for i=1:loop
    fprintf('%d/%d\n',i,loop);
    if isempty(nums)
        [seq(i,:),~,~,oc_map]=pmap(ensemble{i,3},LTraj,PosT,0,'animal','rat','speed',ThS,'verbose',0);
    else
        [seq(i,:),~,~,oc_map]=pmap(extractDelay(ensemble{i,3},event,nums,'jitterPre',jitterPre,'jitterPost',jitterPost),LTraj,PosT,0,'animal','rat', ...
                 'speed',ThS,'verbose',0);
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