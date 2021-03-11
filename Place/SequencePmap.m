function [seq,order,oc_map,LTraj,LTrajS,TMpoints]=SequencePmap(ensemble,event,nums,Traj,PosT,ThS,LTraj,varargin)
if nargin==5
    ThS=2.5;
elseif nargin<=6
    [LTraj,LTrajS,TMpoints]=fitLinear(Traj);
elseif nargin>=7
    LTrajS=[];
    TMpoints=[];
end

p = inputParser;
p.addParamValue('outputmode', 'normalized', @ischar);

p.parse(varargin{:});
outputMode = p.Results.outputmode;


loop=size(ensemble,1);
for i=1:loop
    fprintf('%d/%d\n',i,loop);
    [seq(i,:),~,~,oc_map]=pmap(extractDelay(ensemble{i,3},event,nums),LTraj,PosT,0,'animal','rat', ...
                 'speed',ThS,'verbose',0);
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