function [seq,order,oc_map,LTraj,LTrajS]=SequencePmap(ensemble,event,nums,Traj,PosT,ThS,LTraj)
if nargin==5
    ThS=2.5;
end

if nargin<=6
    [LTraj,LTrajS]=fitLinear(Traj);
end

loop=size(ensemble,1);
for i=1:loop
    fprintf('%d/%d\n',i,loop);
    [seq(i,:),~,~,oc_map]=pmap(extractDelay(ensemble{i,3},event,nums),LTraj,PosT,0,'animal','rat', ...
                 'speed',ThS,'verbose',0);
end

[maxSeq,ind]=max(seq,[],2);
[~,ind]=sort(ind);
seq=seq./maxSeq;
order=ind;

return;