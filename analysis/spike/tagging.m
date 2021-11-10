function minds=tagging(units,tagPoint)

loop=size(units,1);
minds=zeros(1,loop);
for i=1:loop
  minds(i)=tagIndex(units{i},tagPoint);
end

return;
%%%%%%%%%%%%%%%%%%%%%
function mind=tagIndex(unit,tagPoint)

kHz=25;
dispRangePre=kHz*200;%ms before stimOn
dispRangePost=kHz*200;%ms after stimOff
stimDuration=kHz*50;%ms
%stimOff=dispRangePre+kHz*20;
binWidth=1;%ms
step=binWidth/(1/kHz);

triggerPoint=tagPoint;

raster=[];
Ind=[];
cnt=1;

tmps=unit;
tmps=double(tmps);
for i=1:length(triggerPoint)
  raster=[raster tmps-triggerPoint(i)];
  Ind=[Ind ones(size(tmps))*cnt];
  cnt=cnt+1;
end

c=hist(raster,-dispRangePre-step:step:(stimDuration+dispRangePost)+step);
c([1 end])=[];
if sum(c) > 100
  c=(c./(cnt-1))./(binWidth/1000);
  mind=mean(c([dispRangePre/kHz+1:dispRangePre/kHz+stimDuration/kHz]))/mean(c([1:dispRangePre/kHz end-dispRangePost/kHz:end]));
else
  mind=0;
end


return;

