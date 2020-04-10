function [PosT,Traj]=extractPos(filename,Traj)
PosId=8;

[path,name,ext]=fileparts(filename);
Efilename = fullfile(path,[name '/' name 'Event.mat']);

load(Efilename,'event');
Th=2.5;

[PosT]=getTimes(event(PosId,:),Th)';
if size(Traj,1)<size(PosT,1)
    PosT=PosT(1:size(Traj,1),:);
elseif size(Traj,1)>size(PosT,1)
    Traj=Traj(1:size(PosT,1),:);
end


return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%get times
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PosT=getTimes(x,Th)

P=find(x>Th);
P=[0 P];

tmp=find((diff(P)~=1));
tmp=unique([tmp tmp+1]);
tmp(1)=[];

PosT=P(tmp);

%remove oversampled data(interval <5msec)
ro=find(diff(PosT)<125);%125
ro=unique([ro ro+1]);

if size(ro,2)>1
  fprintf('removing oversampled data:%d...\n',size(ro,2));
end
PosT(ro)=[];

if x(1) > Th
    PosT=PosT(2:2:end);
else
    PosT=PosT(1:2:end);
end

return;