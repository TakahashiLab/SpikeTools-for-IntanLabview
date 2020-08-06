function [seq,order]=SequencePmap(ensemble,Traj,PosT,LTraj)
if nargin==3
    [LTraj]=fitLinear(Traj);
end

loop=size(ensemble,1);
for i=1:loop
    fprintf('%d/%d\n',i,loop);
    seq(i,:)=pmap(ensemble{i,3},LTraj,PosT,0,'animal','rat', ...
                 'verbose',0);
end

[maxSeq,ind]=max(seq,[],2);
[~,ind]=sort(ind);
seq=seq(ind,:)./maxSeq(ind);
order=ind;

return;