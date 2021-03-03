%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%get times
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [riseT,fallT]=getTimes(x,Th)

P=find(x>Th);

%if event width==1
P(x(P+1)<Th)=[];

P=[0 P length(x)];


tmp=find((diff(P)~=1));
tmp=unique([tmp tmp+1]);
tmp(1)=[];

PosT=P(tmp);

%remove oversampled data(interval <5msec)
ro=find(diff(PosT)<125);%125:5msec 25000:1sec
ro=unique([ro ro+1]);

 
if size(ro,2)>1
  fprintf('removing oversampled data:%d...\n',size(ro,2));
end

if x(1) > Th
    riseT=PosT(2:2:end-1);
    fallT=PosT(3:2:end);
else
    riseT=PosT(1:2:end-1);
    fallT=PosT(2:2:end);
end

return;

