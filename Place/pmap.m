%Pmap=pmap(spikes,Traj,PosT,0,'animal','rat')
function [rate_map,spatial_scale,SpksShuffle,oc_map,rate_mapB,FrIndex]=pmap(Spks,Traj, msT,fstart,varargin)
verbose=0;
spON=0;
corrFlag=0;
year=2018;
shuffleN=100;
shuffleLag=20000;%20sec
SpksShuffle=[];
BinWidthCm=2.5;%2.5cm
Linear=0;
clockWise=1;
xl=3;
yl=4;

kHz=31.25;
ThS=2;%100cm/s

if size(Traj,1)==1 | size(Traj,2)==1
    Len=max(Traj);
    Linear=1;
    xl=1;
    yl=1;
else
    yLen=max(Traj(:,2))-min(Traj(:,2));
    xLen=max(Traj(:,1))-min(Traj(:,1));
    Len=yLen;
    if yLen < xLen
        Len=xLen;
    end
%     Len=min([xLen yLen])
end


for i=1:2:(length(varargin)-1)
    % change the value of parameter
  switch lower(varargin{i})
    case 'animal',
      if strcmp(varargin{i+1},'bird')
          cmPerPixel=120/Len;%120cm circle;
          Bin=BinWidthCm/cmPerPixel;
          spON=0;

      elseif strcmp(varargin{i+1},'sgbird')
          cmPerPixel=120/Len;%120cm circle;
          Bin=BinWidthCm/cmPerPixel;
          kHz=30;
          spON=0;         

      elseif strcmp(varargin{i+1},'birdsmall')
          cmPerPixel=80/Len;%120cm circle;
          Bin=BinWidthCm/cmPerPixel;
          spON=0;

      elseif strcmp(varargin{i+1},'birdsmalls')
          cmPerPixel=80/Len;%120cm circle;
          Bin=BinWidthCm/cmPerPixel;
          spON=1; 
          ThS=2.5;

      elseif strcmp(varargin{i+1},'birds')
          cmPerPixel=120/Len;%120cm circle;
          Bin=BinWidthCm/cmPerPixel;
          spON=1;
          ThS=2.5;

      elseif strcmp(varargin{i+1},'sgbirds')
          cmPerPixel=120/Len;%120cm circle;
          Bin=BinWidthCm/cmPerPixel;
          kHz=30;
          spON=1;
          ThS=2.5;

      elseif strcmp(varargin{i+1},'sg2birds')
          cmPerPixel=120/Len;%120cm circle;
          Bin=BinWidthCm/cmPerPixel;
          kHz=20;
          spON=1;
          ThS=2.5;

      elseif strcmp(varargin{i+1},'rat')
          if Linear
              cmPerPixel=400/Len;%
          else
              cmPerPixel=160/Len;%50
          end
          Bin=BinWidthCm/cmPerPixel;
          kHz=25;
          spON=1;
          ThS=2.5;        
          msT=msT/kHz;
      end
      shuffle=0;
    case 'shuffle',
      shuffle=1;
    case 'speed',
      ThS=varargin{i+1};
    case 'corr',
      corrFlag=varargin{i+1};
    case 'shufflen',
      shuffleN=varargin{i+1};
    case 'small'
      small=varargin{i+1};
      if small
          cmPerPixel=80/Len;%120cm circle;
      end
    case 'year'
      year=varargin{i+1};
      if year==2019 
          if small
              cmPerPixel=80/Len;%120cm circle;
          else
              cmPerPixel=115/Len;%115cm circle;
          end
      end
    case 'verbose',
      verbose=varargin{i+1};
  end
end

spatial_scale=cmPerPixel;

FPS=floor(1/(median(diff(msT))/1000));
fs_video=FPS;
msFPS=floor(1/FPS*1000);

Spks=ceil(Spks/kHz)+fstart;%msec

if Linear
    dTraj=diff(Traj);
    ind=find(abs(dTraj)>max(Traj)/2);
    dTraj(ind)=dTraj(ind)-sign(dTraj(ind))*max(Traj);
    movement=sqrt(sum(dTraj.^2,2))*cmPerPixel;% 
                                              
else
    movement=sqrt(sum(diff(Traj(:,5:6)).^2,2))*cmPerPixel;% 1:2 -> 5:6
end

speed=movement./diff(msT);
speed=smooth(speed,FPS,'moving')*1000;

%analyzed if the speed is more than 5cm/s 
if spON
  Good= find(speed >= ThS );
  Traj=Traj(Good,:);
  msT=msT(Good);
  speed=speed(Good);
end

if Linear
    x=Traj;
    y=ones(size(Traj));
else
    x=ceil(Traj(:,xl));
    y=ceil(Traj(:,yl));
end


StartTraj=msT(1);
EndTraj=msT(end);
Spks=Spks(find(Spks > StartTraj & Spks < EndTraj));


if isempty(Spks)
    binside=2.5;
    if Linear
        dim=binside;
        xdim = min(x):spatial_scale^-1*dim:max(x); %edges of x and y
        dummy = histcounts(x, xdim);
    else
        dummy = Occupancy(x,y, spatial_scale,binside,fs_video);
    end
    
    rate_map=zeros(size(dummy));
    spatial_scale=[];
    SpksShuffle=[];
    oc_map=ones(size(dummy));
    rate_mapB=[];
    FrIndex=[];
    return;
end

%Firing rate contrast Index
windowSize=2000;%2s
Th=nanmedian(speed);
FrHigh=[];
FrLow=[];
for k=StartTraj:windowSize:(EndTraj-windowSize)
    fr=sum(Spks > k & Spks < k+windowSize);
    instSpeed=nanmedian(speed(msT > k & msT < k+windowSize));
    
    if instSpeed>Th
        FrHigh=[FrHigh fr];
    else
        FrLow=[FrLow fr];
    end
end

FrHigh=sum(FrHigh)/(length(FrHigh)*(windowSize/1000));
FrLow=sum(FrLow)/(length(FrLow)*(windowSize/1000));

FrIndex(1,1)=(FrHigh-FrLow)/(FrHigh+FrLow);
FrIndex(1,2)=FrHigh;
FrIndex(1,3)=FrLow;


spk_x=[];
spk_y=[];

fieldShuffleFlag=0;

if shuffle
    
    beginSpks=msT(1);
    endSpks=msT(end);
    entireLength=ceil(endSpks-beginSpks);
    rp=randperm(entireLength-shuffleLag*2+1);
    rp=rp(1:shuffleN);
    orgRp=shuffleLag:(entireLength-shuffleLag);
    rp=orgRp(rp);
    lenOrgRp=length(orgRp);
    lenSpks=length(Spks);
    SpksShuffle=zeros(shuffleN,lenSpks);
    for i=1:shuffleN
        spkss=Spks+rp(i);
        topSpks=find(spkss>=endSpks);
        spkss(topSpks)=spkss(topSpks)-endSpks+beginSpks;
        SpksShuffle(i,:)=sort(spkss);
    end
    
    if ~fieldShuffleFlag
        for i=1:shuffleN
            Spks=SpksShuffle(i,:);
            spk_x=[];
            spk_y=[];
            for j=1:(size(Traj,1)-1)
                SpkCnt=sum(Spks >= msT(j) & Spks < msT(j)+msFPS);
                for k=1:SpkCnt
                    spk_x=[spk_x; Traj(j,xl)];
                    spk_y=[spk_y; Traj(j,yl)];
                end
            end

            
            if Linear
                [rate_map(i,:,:)] = rateLinearMap(spk_x,x,spatial_scale,fs_video,1);
            else
                [rate_map(i,:,:)] = ratemap(spk_x,spk_y,x,y, ...
                                            spatial_scale,fs_video);
            end
        end
        oc_map=[];
    else%fieldShuffle
    
        SpkCnt=[];
        spk_x=[];
        spk_y=[];
        for j=1:(size(Traj,1)-1)
            spkcnt=sum(Spks >= msT(j) & Spks < msT(j)+msFPS);
            if spkcnt
                for k=1:spkcnt
                    spk_x=[spk_x; Traj(j,xl)];
                    spk_y=[spk_y; Traj(j,yl)];
                end
            end
        end

        [rm,~,~,oc_map] = ratemap(spk_x,spk_y,x,y,spatial_scale,fs_video,1);
        

        for i=1:shuffleN
            [rate_map(i,:,:)] = fieldShuffleEE(rm,oc_map);
        end
        
        %for Border cell
        for i=1:shuffleN
            Spks=SpksShuffle(i,:);
            spk_x=[];
            spk_y=[];
            for j=1:(size(Traj,1)-1)
                SpkCnt=sum(Spks >= msT(j) & Spks < msT(j)+msFPS);
                for k=1:SpkCnt
                    spk_x=[spk_x; Traj(j,xl)];
                    spk_y=[spk_y; Traj(j,yl)];
                end
            end
            [rate_mapB(i,:,:)] = ratemap(spk_x,spk_y,x,y,spatial_scale,fs_video);
        end
        
    end
else
    SpkCnt=[];
    spk_x=[];
    spk_y=[];
    spk_x_disp=[];
    spk_y_disp=[];
    jitter=[-5:1:5];
    
    for j=1:(size(Traj,1)-1)
        spkcnt=sum(Spks >= msT(j) & Spks < msT(j)+msFPS);

        if spkcnt
            for k=1:spkcnt
                rp=randperm(11);
                jitterS=jitter(rp(1));
                spk_x_disp=[spk_x_disp; Traj(j,xl)+jitterS];
                spk_y_disp=[spk_y_disp; Traj(j,yl)+jitterS];
                spk_x=[spk_x; Traj(j,xl)];
                spk_y=[spk_y; Traj(j,yl)];
            end
        end
    end

    if verbose
        line(Traj(:,xl),Traj(:,yl),'Color','k');
        hold on;
        scatter(spk_x,spk_y,'.r','SizeData',25);
        axis equal off;
        set(gca,'xdir','normal');
        set(gca,'ydir','reverse');
        maxX=max(spk_x);
        maxY=max(spk_y);
    end
  
    if Linear
        [rate_map,~,oc_map] = rateLinearMap(spk_x,x,spatial_scale,fs_video,1);
     
    else
        [rate_map,~,~,oc_map] = ratemap(spk_x,spk_y,x,y, ...
                                        spatial_scale,fs_video,1);
    end
        rate_mapB=[];
end

