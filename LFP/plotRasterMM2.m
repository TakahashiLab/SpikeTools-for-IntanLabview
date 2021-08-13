%verbose figureを出したくない時は0
%load event.mat
%plotRaster(kkOut{4,3},event,[1 2],'binWidth',100);
%%% Meanがevent=[1 2]で作った時と各event別で作った時で数字が合わない
%Nsta,Nen:解析したい周回数のstart,end
%SpkInIndex:TM中の活動
function [SpkInIndex,histR,Up,Do,smhistR]=plotRasterMM2(spks,event,nums,PORK3,Nsta,Nen,varargin)
p = inputParser;
p.addParamValue('binwidth', 100, @isnumeric);
p.addParamValue('samplingrate', 25, @isnumeric);
p.addParamValue('jitter', 0, @isnumeric);
p.addParamValue('smooth', 10, @isnumeric);
p.addParamValue('verbose', 1, @isnumeric);
p.addParamValue('splitNum', 0, @isnumeric);
p.addParamValue('phase', 0, @isnumeric);

p.parse(varargin{:});
binWidth = p.Results.binwidth;
kHz=p.Results.samplingrate;
jitter=p.Results.jitter;
sm = p.Results.smooth;
verbose = p.Results.verbose;
splitNum = p.Results.splitNum;
calcPhase = p.Results.phase;

phaseT=spks{1,4};
spks=double(spks{1,3});


MAX=0;
Th=1.5;%event threshold
ThS=5;%split threshold (sec)
loop=length(nums);
SpkInIndex=[];

for m=nums
    if verbose
          if loop~=1
            subplot(1,loop,m);
          end
    end
    TM=find((diff(event(m,:))>Th));
    TM(find(((event(m,TM+2)))<Th))=[];
    TM_fall=find((diff(event(m,:))<-Th));
    TM_fall(find(((event(m,TM_fall-1)))<Th))=[];
    Up=[TM(find(PORK3(Nsta)<TM & PORK3(Nen)>TM))];
    Do=[TM_fall(find(PORK3(Nsta)<TM_fall & PORK3(Nen)>TM_fall))];
    l=size(Up,2);

%     
%     [Up,Do]=getTimes(event(nums(m),:),Th);
%     
% 
%     if event(nums(m),1)>Th
%         buf=Do;
%         Do=Up;
%         Up=buf;
%     end
%     
%     if Up(1)> Do(1)
%        Up=Up(1:end-1);
%        Do=Do(2:end);
%        l=length(Up);       
%     elseif length(Up)~=length(Do)
%        l=min([length(Up) length(Do)]);
%        Up=Up(1:l);
%        Do=Do(1:l);
%     else
%       l=length(Up);
%     end
    
    
%     if splitNum>0
%        [Up,Do,l]=splitUpDo(Up,Do,splitNum,ThS,kHz);
%     end

    duration=max(Do-Up)+jitter*1000*kHz*2;
    eventLength=floor(duration/(kHz*binWidth));
    if calcPhase
        raster=zeros(360,eventLength);
    else
        raster=zeros(l,eventLength);
    end
    
    if verbose
        hold on;
    end
    for i=1:l
        SpkInTrial=spks-(Up(i)-jitter*1000*kHz);
        SpkInIndex=[SpkInIndex find(SpkInTrial>=0 & SpkInTrial<= ...
                                    duration)];
        SpkInIndexT=find(SpkInTrial>=0 & SpkInTrial<=duration);        
        SpkInTrial=SpkInTrial(SpkInTrial>=0 & SpkInTrial<=duration);
        SpkSeq=floor(SpkInTrial/(kHz*binWidth));
        SpkInIndexT(SpkSeq==0)=[];
        SpkSeq(SpkSeq==0)=[];
        SpkInIndexT(SpkSeq==eventLength+1)=[];
        SpkSeq(SpkSeq==eventLength+1)=[];
        
        if calcPhase
            phaseSeq=phaseT(SpkInIndexT);
            raster=raster+full(sparse(phaseSeq,SpkSeq,1,360,eventLength));
        else
            raster(i,:)=full(sparse(1,SpkSeq,1,1,eventLength));
            spkPoint=find(raster(i,:));
            if ~isempty(spkPoint)
                if verbose
                    plot(spkPoint,i,'k.');
                end
            end
        end
    end
    

    if calcPhase & verbose
        raster=[raster;raster];
        [x,y]=find(raster);
        plot(y,x,'k.');
    end
    
    histR{m}=sum(raster)/l;
    MAX=max([MAX max(histR{m})]);
    
    jump=1000/binWidth;
    if verbose

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
    
    if calcPhase
        axis([0 size(raster,2) 0 size(raster,1)]);
        ylabel('Theta phase (degree)');
    else
        ylabel('Lap #');
    end
    end
end



    for m=nums
        
        if verbose
            if loop~=1
                subplot(1,loop,m);
            end
            hold on;
        end
        
        histRaster=histR{m};
       
        fr=(max(histRaster)/(binWidth/1000));
        histRaster=histRaster/MAX*l;
        if sm
           smhistRaster=smooth(histRaster,sm);
        else
            smhistRaster=histRaster;
        end
%         MAX=max([MAX max(smhistRaster)]);
%         histRaster=histRaster/MAX*l;
%         smhistRaster=smhistRaster/MAX*(l+1);
        [~,Index]=max(smhistRaster);%mieno
        index=Index/jump-jitter;
        Mean=mean(histRaster);%mieno
        
        if verbose 
            if calcPhase
                plot(smhistRaster./max(smhistRaster)*720,'r');
            else
                plot(smhistRaster,'r');
            end
        title(sprintf('Peak:%2.1fHz Index:%2.1f Mean:%2.1fHz',fr,index,Mean),'FontSize',15);%mieno
        end
        smhistR{m}=smhistRaster;
    end

return;
%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Up,Do,l]=splitUpDo(Up,Do,splitNum,ThS,kHz)
% duration=(Do-Up)/(kHz*1000);
% splitP=find(abs(diff(duration))>ThS )+1;
% splitP=[1 splitP];
% %duration
% if length(splitP)==splitNum
%   r=splitP(splitNum):length(Do);
% else
%   r=splitP(splitNum):(splitP(splitNum+1)-1);
% end
% 
% Up=Up(r);
% Do=Do(r);
% l=length(Up);
% 
% return;