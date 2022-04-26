%%%%%%%%%%%%%%%%%
%Traj/Pos=Posistion(n x posture)
%msT/PosT=time of position data
%posture: Index of head position in Traj array
%fstart: file start time (for logger), not necessary if thethered recording 
%maxdist: the maximum distance of arena/maze in a centimeter scale
%khz:sampling rate (intan(takahashilab):25.0, mouselog:31.25)
%spon: speed thresholding
%speed: cm/s
%binwidth: place field binning in a centimeter scale
%grid: kernel size
%verbose: output
%%%%%%%%%%%%%%%%
%%For bird,
%[ratemap,~,~,ocmap]=pmap2(ensemble{i,3},Traj,msT,'fstart',fstart,'animal','any','maxdist',150,'khz',31.25,'spon',1,'speed',5,'binwidth',2.5,'grid',[10 10],'verbose',1);
%%%%%%%%%%%%%%%%
%%For mouse,%
%[ratemap,~,~,ocmap]=pmap2(ensemble{i,3},Traj,msT,'animal','rodent','maxdist',150,'khz',25,'spon',1,'speed',2.5,'binwidth',2.5,'grid',[10
%10]);
%
%
%%for display
%%imagePmap(ratemap,ocmap);


function [rate_map,spatial_scale,SpksShuffle,oc_map,rate_mapB,FrIndex,xdim,ydim]=pmap2(Spks,Traj, msT,varargin)
corrFlag=0;
year=2018;
shuffleN=100;
shuffleLag=20000;%20sec
SpksShuffle=[];
clockWise=1;
xl=3;
yl=4;

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

p = inputParser;
p.addParamValue('fstart', 0, @isnumeric);
p.addParamValue('linear', 0, @isnumeric);
p.addParamValue('maxdist', 400, @isnumeric);
p.addParamValue('animal', 'any', @ischar);
p.addParamValue('khz', 31.25, @isnumeric);
p.addParamValue('spon', 0, @isnumeric);
p.addParamValue('speed', 2.5, @isnumeric);
p.addParamValue('shuffle', 0, @isnumeric);
p.addParamValue('shufflen', 0, @isnumeric);
p.addParamValue('verbose', 0, @isnumeric);
p.addParamValue('posture', [3 4 5 6], @isvector);
p.addParamValue('binwidth', 2.5, @isnumeric);
p.addParamValue('std', 3, @isnumeric);
p.addParamValue('xdim', [], @isvector);
p.addParamValue('ydim', [], @isvector);
p.addParamValue('grid', [], @isvector);
p.addParamValue('spatialscale', -1, @isnumeric);


p.parse(varargin{:});
fstart = p.Results.fstart;
Linear = p.Results.linear;
maxdist=p.Results.maxdist;
animal=p.Results.animal;
kHz=p.Results.khz;
spON = p.Results.spon;
ThS=p.Results.speed;
shuffle=p.Results.shuffle;
shuffleN=p.Results.shufflen;
verbose=p.Results.verbose;
posture = p.Results.posture;
BinWidthCm=p.Results.binwidth;
binside=p.Results.std;
spatial_scale=p.Results.spatialscale;


xdim = p.Results.xdim;
ydim = p.Results.ydim;
GRID = p.Results.grid;

cmPerPixel=maxdist/Len;%

%%%%%%%%%%%%%%
%%%For bird, large arena 
%%'animal','bird','maxdist',120
%%
%%if spikegadgets
%%%'kHz',30
%%if small arena
%%'maxdist',80
%%%%%%%%%%%
%%%%%%%%%%%%%%
%%%For rat, 
%%%if Linear, then 'maxdist',400
%%%else 'maxdist',160


if size(Traj,1)~=size(msT,1)
    msFPS=median(diff(msT));
    seq=1:length(Traj);
    seq=seq-1;
    msT=seq*msFPS+msT(1);
    msT=msT';
end


switch lower(animal)
  case 'any'
    cmPerPixel=maxdist/Len;%120cm circle;
  case 'rodent'
    cmPerPixel=maxdist/Len;%
    msT=msT/kHz;
end

if spatial_scale<0
    spatial_scale=cmPerPixel/BinWidthCm;
end

FPS=floor(1/(median(diff(msT))/1000));
fs_video=FPS;
msFPS=floor(1/FPS*1000);

Spks=ceil(Spks/kHz)+fstart;%msec

Traj=Traj';
x=Traj(posture(1),:);
y=Traj(posture(2),:);
x2=Traj(posture(3),:);
y2=Traj(posture(4),:);
Traj=Traj';
headdir=atan2d(x-x2,y-y2)';

%Good=find( headdir>0 & headdir <90);
%Good=find( headdir>90);
%Good=find( headdir<0 & headdir>-90);
Good=find( headdir<-90);

%Traj=Traj(Good,:);
%msT=msT(Good);

if Linear
    dTraj=diff(Traj);
    ind=find(abs(dTraj)>max(Traj)/2);
    dTraj(ind)=dTraj(ind)-sign(dTraj(ind))*max(Traj);
    movement=sqrt(sum(dTraj.^2,2))*cmPerPixel;% 
                                              
else
    movement=sqrt(sum(diff(Traj(:,posture(3:4))).^2,2))*cmPerPixel;% 1:2 -> 5:6
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
                [rate_map(i,:,:)] = rateLinearMap(spk_x,x,spatial_scale,fs_video,binside);
            else
                [rate_map(i,:,:)] = ratemap(spk_x,spk_y,x,y, ...
                                            spatial_scale,fs_video,binside);
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

        [rm,~,~,oc_map] = ratemap(spk_x,spk_y,x,y,spatial_scale,fs_video,binside);
        

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
            [rate_mapB(i,:,:)] = ratemap(spk_x,spk_y,x,y,spatial_scale,fs_video,binside);
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
        [rate_map,xdim,oc_map] = rateLinearMap(spk_x,x,spatial_scale,fs_video,binside);
     
    else
        if ~isempty(xdim) & ~isempty(GRID)
            [rate_map,xdim,ydim,oc_map] = ratemap(spk_x,spk_y,x,y, ...
                                                  spatial_scale, ...
                                                  fs_video,binside,GRID,xdim,ydim);
        elseif  ~isempty(GRID)
            [rate_map,xdim,ydim,oc_map] = ratemap(spk_x,spk_y,x,y, ...
                                                  spatial_scale, ...
                                                  fs_video,binside,GRID);        
        else
            [rate_map,xdim,ydim,oc_map] = ratemap(spk_x,spk_y,x,y, ...
                                                  spatial_scale, ...
                                                  fs_video,binside);            
        end
        
    end
        rate_mapB=[];
end

