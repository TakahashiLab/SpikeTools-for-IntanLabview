%
%load event.mat
%plotRaster(kkOut{4,3},event,[1 2],'binWidth',100);
%%%
function histR=plotRasterL(spks,event,num,varargin)
p = inputParser;
p.addParamValue('binwidth', 100, @isnumeric);
p.addParamValue('samplingrate', 25, @isnumeric);
p.addParamValue('jitter', 0, @isnumeric);
p.addParamValue('smooth', 0, @isnumeric);
p.addParamValue('verbose', 1, @isnumeric);

p.parse(varargin{:});
binWidth = p.Results.binwidth;
kHz=p.Results.samplingrate;
jitter=p.Results.jitter;
sm = p.Results.smooth;
verbose = p.Results.verbose;

spks=double(spks);

MAX=0;
Th=2.5;%event threshold

nums=[num num];
loop=2;%left,right direction

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
    end
    l=ceil(l/2);
    eventNum=l;
    
    if m==1
      range=1:2:l;
    else
      range=2:2:l;
    end
    
    duration=max(Do-Up)+jitter*1000*kHz*2;
    eventLength=floor(duration/(kHz*binWidth));
    raster=zeros(l,eventLength);
    c=1;
    for i=range
        SpkInTrial=spks-(Up(i)-jitter*1000*kHz);
        SpkInTrial=SpkInTrial(SpkInTrial>=0 & SpkInTrial<=duration);
        SpkSeq=floor(SpkInTrial/(kHz*binWidth));
        SpkSeq(SpkSeq==0)=1;
        raster(c,:)=full(sparse(1,SpkSeq,1,1,eventLength));
        
        spkPoint=find(raster(c,:));
        if ~isempty(spkPoint)
            if verbose
               hold on;
                plot(spkPoint,c,'k.');
            end
        end
        c=c+1;
    end

    histR{m}=sum(raster);
    MAX=max([MAX max(histR{m})]);
    
    if verbose
    jump=1000/binWidth;
    xticks=[0:jump:duration/(kHz*binWidth)];
    xlabels=floor(xticks/(jump));
    
    if jitter
        xlabels=xlabels-jitter;
        plot([jitter*jump jitter*jump],[0 l+1],'g')
        plot([duration/(kHz*binWidth)-jitter*jump duration/(kHz*binWidth)-jitter*jump],[0 l+1],'g')
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
            plot(smooth(histRaster,10),'r');
        else
            plot(histRaster,'r');
        end
        title(sprintf('Peak:%2.1fHz',fr));
    end
end
return;