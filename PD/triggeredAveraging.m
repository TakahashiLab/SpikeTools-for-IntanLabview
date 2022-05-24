function ta=triggeredAveraging(x,tagPoint)

kHz=25;
duration=10*kHz;%ms;
Th=200;

loop=size(tagPoint,2);
ta=zeros(duration*2+1,loop);
for i=1:loop
  pre=tagPoint(i)-duration;
  post=tagPoint(i)+duration;
  ta(:,i)=x(pre:post);
  ta(ta(:,i)>-Th  & ta(:,i)<Th,i)=0;
end


return;
