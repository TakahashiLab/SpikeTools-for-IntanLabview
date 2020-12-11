function [PosT,Traj]=extractPos(filename,Traj,PosId,ver)
if nargin==2
    PosId=8;
    ver=1;%new
elseif nargin==3
    ver=1;
end

[path,name,ext]=fileparts(filename);
Efilename = fullfile(path,[name '/' name 'Event.mat']);

load(Efilename,'event');
Th=2.5;

[PosT]=getTimes(event(PosId,:),Th,ver)';
if length(Traj)<length(PosT)
    PosT=PosT(1:length(Traj),:);
elseif length(Traj)>length(PosT)
    Traj=Traj(1:length(PosT),:);
end


return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%get times
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PosT=getTimes(x,Th,ver)

P=find(x>Th);
P=[0 P];

tmp=find((diff(P)~=1));
tmp=unique([tmp tmp+1]);


tmp(1)=[];

PosT=P(tmp);
size(PosT)

%remove oversampled data(interval <5msec)
ro=find(diff(PosT)<125);%125
ro=unique([ro ro+1]);

if size(ro,2)>1
  fprintf('removing oversampled data:%d...\n',size(ro,2));
end
size(PosT)
PosT(ro)=[];

if ver
    size(PosT)
    if x(1) > Th
        PosT=PosT(2:2:end);
    else
        PosT=PosT(1:2:end);
    end
    size(PosT)
else
    dp=diff(PosT);
    PosT(find(dp<=463)+1)=[];
end

return;