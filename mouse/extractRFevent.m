function segPara=extractRFevent(TrialT,Pos,PosT)
kHz=25;
loop=size(TrialT,2)+1;

pre=[1 TrialT+1];
post=[TrialT PosT(end)];

segments=[pre' post'];

segPara=cell(loop,2);
%for i=1
for i=1:loop
  fprintf('seg #%d\n',i);
  seg=find(PosT>segments(i,1) & PosT<segments(i,2));
  pos=Pos(seg,1:2);  
  post=PosT(seg);  
  [segPara{i,1},segPara{i,2}]=runSpeedMaxDist(pos,post,kHz);%running speed /sec
end

return;
%%%%%%%%%%%%%%%%%%%%%%%%
function [speed,maxDist]=runSpeedMaxDist(pos,post,kHz)

halfsec=kHz*1000/2;

dist=sqrt(sum(diff(pos).^2,2));
dt=diff(post);
%endP=size(pos,1);
endP=length(post)-1;
loop=length(post);
speed=[];
maxDist=[];

if post(2)-post(1)> 1000
    post(1)=[];
end

for i=1:60:loop
  seg=find(post>post(i)-halfsec & post<post(i)+halfsec);
  if seg(end)>endP
    seg(end)=[];
  end
%  seg(find(sum(isnan(pos(seg,:)),2)))=[];

  dpos=pos([seg(1) seg(end)],:);
  %speed=[speed sqrt(sum(diff(dpos).^2,2))];
  speed=[speed nanmean(dist(seg)./dt(seg))];
  p=pos(seg)-pos(i);
  maxDist=[maxDist nanmax(sqrt(sum(diff(p).^2,2)))];
end


return;

