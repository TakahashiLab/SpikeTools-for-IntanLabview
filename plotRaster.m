%
%load event.mat
%plotRaster(kkOut{4,3},event,[1 2],'binWidth',100);
%%%
function [histR,Up,Do,I,M]=plotRaster(spks,event,nums,varargin)
p = inputParser;
p.addParamValue('binwidth', 100, @isnumeric);
p.addParamValue('samplingrate', 25, @isnumeric);
p.addParamValue('jitter', 0, @isnumeric);
p.addParamValue('smooth', 0, @isnumeric);
p.addParamValue('verbose', 1, @isnumeric);
p.addParamValue('splitNum', 0, @isnumeric);

p.parse(varargin{:});
binWidth = p.Results.binwidth;
kHz=p.Results.samplingrate;
jitter=p.Results.jitter;
sm = p.Results.smooth;
verbose = p.Results.verbose;
splitNum = p.Results.splitNum;

spks=double(spks);


MAX=0;
Th=3;%event threshold
ThS=5;%split threshold (sec)
loop=length(nums);
for m=1:loop
    if verbose
        subplot(1,loop,m);
    end
    
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
       l=min([length(Up) length(Do)]);
       Up=Up(1:l);
       Do=Do(1:l);
    else
      l=length(Up);
    end
    
    
    if splitNum>0
       [Up,Do,l]=splitUpDo(Up,Do,splitNum,ThS,kHz);
    end
    
   
      
    eventNum=l;
    duration=max(Do-Up)+jitter*1000*kHz*2;
    eventLength=floor(duration/(kHz*binWidth));
    raster=zeros(l,eventLength);
    if verbose
        hold on;
    end
    for i=1:l
        SpkInTrial=spks-(Up(i)-jitter*1000*kHz);
        SpkInTrial=SpkInTrial(SpkInTrial>=0 & SpkInTrial<=duration);
        SpkSeq=floor(SpkInTrial/(kHz*binWidth));
        SpkSeq(SpkSeq==0)=[];
        SpkSeq(SpkSeq==eventLength+1)=[];
        raster(i,:)=full(sparse(1,SpkSeq,1,1,eventLength));
        
        spkPoint=find(raster(i,:));
        if ~isempty(spkPoint)
            if verbose
                plot(spkPoint,i,'k.');
            end
        end
    end

    histR{m}=sum(raster)/l;
    MAX=max([MAX max(histR{m})]);
    
    if verbose
    jump=1000/binWidth;
    xticks=[0:jump:duration/(kHz*binWidth)];
    xlabels=floor(xticks/(jump));
    
    if jitter
        xlabels=xlabels-jitter;
        if jitter<0
          pjitter=0;
        else
          pjitter=jitter;
        end
          plot([pjitter*jump pjitter*jump],[0 l+1],'g')
          plot([duration/(kHz*binWidth)-pjitter*jump duration/(kHz*binWidth)-pjitter*jump],[0 l+1],'g')
        
    end

    set(gca,'xtick',xticks,'xticklabels',xlabels);
    xlabel('Delay time (s)');
    ylabel('Lap #');
    end
end

if verbose
    for m=1:loop
        subplot(1,loop,m);
        hold on;
        histRaster=histR{m};
       
        fr=(max(histRaster)/(binWidth/1000))/l;
        histRaster=histRaster/MAX*l;
        if sm
           smhistRaster=smooth(histRaster,10);
        else
            smhistRaster=histRaster;
        end
        [~,Index]=max(smhistRaster);%mieno
        index=Index/jump-jitter;
        Mean=mean(histRaster);%mieno
        plot(smhistRaster,'r');
        title(sprintf('Peak:%2.1fHz Index:%2.1f Mean:%2.1f',fr,index,Mean));%mieno
        I(1,m)=index;
        M(1,m)=Mean;
    end
end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%
function [Up,Do,l]=splitUpDo(Up,Do,splitNum,ThS,kHz)
duration=(Do-Up)/(kHz*1000);
splitP=find(abs(diff(duration))>ThS )+1;
splitP=[1 splitP];
%duration
if length(splitP)==splitNum
  r=splitP(splitNum):length(Do);
else
  r=splitP(splitNum):(splitP(splitNum+1)-1);
end

Up=Up(r);
Do=Do(r);
l=length(Up);

return;