%
%load event.mat
%plotRaster(kkOut{4,3},event,[1 2],'binWidth',100);
%%%
function [spks,Up,Do]=extractDelay(spks,event,nums,varargin)
p = inputParser;
p.addParamValue('binwidth', 100, @isnumeric);
p.addParamValue('samplingrate', 25, @isnumeric);
p.addParamValue('jitterpre', 0, @isnumeric);
p.addParamValue('jitterpost', 0, @isnumeric);
p.addParamValue('smooth', 0, @isnumeric);

p.parse(varargin{:});
binWidth = p.Results.binwidth;
kHz=p.Results.samplingrate;
jitterPre=p.Results.jitterpre;
jitterPost=p.Results.jitterpost;
sm = p.Results.smooth;

jitterPre=jitterPre*1000*kHz;
jitterPost=jitterPost*1000*kHz;

spks=double(spks);
orgSpks=spks;
MAX=0;
Th=3;%event threshold

delInd=[];
loop=length(nums);
for m=1:loop
    [Up,Do]=getTimes(event(nums(m),:),Th);
    
    if event(nums(m),1)>Th
        buf=Do;
        Do=Up;
        Up=buf;
    end
    
    if Up(1)> Do(1)
       Up=Up(1:end-1);
       Do=Do(2:end);
       l=length(Up);
    elseif length(Up)~=length(Do)
       l=max([length(Up) length(Do)]);
       
       %Up=Up(1:l);
       %Do=Do(1:l);
       Do=[ Do realmax];
    else
        l=length(Up);
    end
    
    eventNum=l;
    for i=l:-1:1
         delInd=[delInd find(orgSpks>=Up(i)-jitterPre & orgSpks<=Do(i)+jitterPost)];
    end
    
end

spks(delInd)=[];
return;