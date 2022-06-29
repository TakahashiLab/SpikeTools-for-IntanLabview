%
%load event.mat
%plotRaster(kkOut{4,3},event,[1 2],'binWidth',100);
%%%
function [histR]=plotTimePhase(spks,Event,phase,varargin)
p = inputParser;
p.addParamValue('binwidth', 100, @isnumeric);
p.addParamValue('samplingrate', 25, @isnumeric);
p.addParamValue('jitterpre', 0, @isnumeric);
p.addParamValue('jitterpost', 0, @isnumeric);
p.addParamValue('smooth', 0, @isnumeric);
p.addParamValue('verbose', 1, @isnumeric);
p.addParamValue('splitNum', 0, @isnumeric);

p.parse(varargin{:});
binWidth = p.Results.binwidth;
kHz=p.Results.samplingrate;
pre=p.Results.jitterpre;
post=p.Results.jitterpost;
sm = p.Results.smooth;
verbose = p.Results.verbose;
splitNum = p.Results.splitNum;

spks=double(spks);
phs=phase(ceil(spks/kHz));
%duration=jitter*1000*kHz*2;
pre=pre*1000*kHz;
post=post*1000*kHz;
duration=pre+post;
eventLength=floor(duration/(kHz*binWidth));
l=length(Event);
raster=zeros(l,eventLength);

if verbose
    hold on;
end

raster=zeros(360,eventLength);
for i=1:l
    SpkInTrial=spks-(Event(i)-pre);
    Id=find(SpkInTrial>=0 & SpkInTrial<=duration);
    SpkInTrial=SpkInTrial(Id);
    SpkPhase=phs';
    SpkPhase=SpkPhase(Id);

    SpkSeq=floor(SpkInTrial/(kHz*binWidth));
    Id=find(SpkSeq==0);
    SpkSeq(Id)=[];
    SpkPhase(Id)=[];

    Id=find(SpkSeq==eventLength+1);
    SpkSeq(Id)=[];
    SpkPhase(Id)=[];
    
    spkPoint=full(sparse(SpkPhase,SpkSeq,1,360,eventLength));
    raster=raster+spkPoint;

    [y,x]=find(spkPoint);
    if ~isempty(spkPoint)
        if verbose
            plot(x,y,'k.');
            plot(x,y+360,'k.');
        end
    end
    

end

if verbose

    histRaster=sum(raster)/l;    
    fr=(max(histRaster)/(binWidth/1000));    
    jump=1000/binWidth;
        xticks=[0:jump:duration/(kHz*binWidth)];
        xlabels=floor(xticks/(jump));

        [~,Index]=max(histRaster);%mieno
        index=Index/jump-pre;
    
        if 1
            pjitter=pre/binWidth/kHz;
            xlabels=xlabels-pjitter/10;
            plot([pjitter pjitter],[0 360*2],'g')
            plot([duration/(kHz*binWidth)-pjitter duration/(kHz*binWidth)-pjitter],[0 360*2],'g')
        
        end

        set(gca,'xtick',xticks,'xticklabels',xlabels);
        xlabel('Delay time (s)');
        ylabel('Theta phase');
    
    title(sprintf('Peak:%2.1fHz Index:%2.1f',fr,index));%
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