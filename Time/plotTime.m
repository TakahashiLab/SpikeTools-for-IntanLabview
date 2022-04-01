%
%load event.mat
%plotRaster(kkOut{4,3},event,[1 2],'binWidth',100);
%%%
function [histR]=plotTime(spks,Event,varargin)
p = inputParser;
p.addParamValue('binwidth', 100, @isnumeric);
p.addParamValue('samplingrate', 25, @isnumeric);
p.addParamValue('pre', 0, @isnumeric);
p.addParamValue('post', 0, @isnumeric);
p.addParamValue('smooth', 0, @isnumeric);
p.addParamValue('verbose', 1, @isnumeric);
p.addParamValue('splitNum', 0, @isnumeric);

p.parse(varargin{:});
binWidth = p.Results.binwidth;
kHz=p.Results.samplingrate;
pre=p.Results.pre;
post=p.Results.post;
sm = p.Results.smooth;
verbose = p.Results.verbose;
splitNum = p.Results.splitNum;

spks=double(spks);
%duration=jitter*1000*kHz*2;
pre=pre*1000*kHz;
post=post*1000*kHz;
duration=pre+post;
eventLength=floor(duration/(kHz*binWidth));
l=length(Event);
raster=zeros(l,eventLength);

hold on;
for i=1:l
    SpkInTrial=spks-(Event(i)-pre);
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
    

    histR=sum(raster)/l;
    MAX=max(histR);
    

    if verbose
        jump=1000/binWidth;
        xticks=[0:jump:duration/(kHz*binWidth)];
        xlabels=floor(xticks/(jump));
    
        if 1
            pjitter=pre/binWidth/kHz;
            xlabels=xlabels-pjitter/10;
            plot([pjitter pjitter],[0 l+1],'g')
            plot([duration/(kHz*binWidth)-pjitter duration/(kHz*binWidth)-pjitter],[0 l+1],'g')
        
        end

        set(gca,'xtick',xticks,'xticklabels',xlabels);
        xlabel('Delay time (s)');
        ylabel('Lap #');
    
    end
end



if verbose

        histRaster=histR;
        fr=(max(histRaster)/(binWidth/1000));

        histRaster=histRaster/MAX*l;
        if sm
           smhistRaster=smooth(histRaster,10);
        else
            smhistRaster=histRaster;
        end
        [~,Index]=max(smhistRaster);%mieno
        index=Index/jump-pre;
        %        [Mean]=mean(histRaster);%mieno

        plot(smhistRaster,'r');
        title(sprintf('Peak:%2.1fHz Index:%2.1f',fr,index));%
        %        title(sprintf('Peak:%2.1fHz Index:%2.1f Mean:%2.1f',fr,index,Mean));%mieno
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