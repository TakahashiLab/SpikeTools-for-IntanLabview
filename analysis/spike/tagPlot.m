function tagPlot(unit,tagPoint)

kHz=25;
dispRangePre=kHz*200;%ms before stimOn
dispRangePost=kHz*200;%ms after stimOff
stimDuration=kHz*50;%ms
%stimOff=dispRangePre+kHz*20;
binWidth=1;%ms
step=binWidth/(1/kHz);

triggerPoint=tagPoint;

clf;
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

stimAfter=dispRangePre+stimDuration;

subplot(2,1,1);
plot(raster,Ind,'.');
axis([-dispRangePre stimDuration+dispRangePost 0 cnt]);
tics=[-dispRangePre 0 stimDuration stimDuration+dispRangePost];
hold on;
plot([0 0],[0 cnt],'r');
plot([stimDuration stimDuration],[0 cnt],'g');
set(gca,'xtick',tics,'xticklabel',tics/kHz);

subplot(2,1,2);
raster(raster<-dispRangePre)=[];
raster(raster>(stimDuration+dispRangePost))=[];

c=hist(raster,-dispRangePre:step:(stimDuration+dispRangePost));
c=(c./(cnt-1))./(binWidth/1000);
bar(c);

bf=mean(c(1:dispRangePre/step));
sf=mean(c(dispRangePre/step:stimAfter/step));

fprintf('base fr=%f, stim fr=%f\n',bf,sf);

maxC=max(c);
if maxC==0
  maxC=1;
end
axis([0 (dispRangePre+stimDuration+dispRangePost)/step 0 maxC*1.2]);
hold on;

plot([dispRangePre/step dispRangePre/step],[0 max(c)*1.2],'r');
plot([stimAfter/step stimAfter/step],[0 max(c)*1.2],'g');

tics=[0 dispRangePre/step stimAfter/step];
ticlabel=[-dispRangePre/step 0 stimDuration/step];
set(gca,'xtick',tics,'xticklabel',ticlabel);

return;
kHz=25;
window=500;%msec
win=window*kHz;

loop=size(tagPoint,2);
unit=unit(unit>tagPoint(1)-win & unit<tagPoint(end)+win);
Spks=unit;
trials=unit;
for i=1:loop
  spk=unit-tagPoint(i);
  ind=find(spk>-win & spk<win);
  Spks(ind)=unit(ind)-tagPoint(i);
  trials(ind)=i;
end


plot(Spks,trials,'.');
return;