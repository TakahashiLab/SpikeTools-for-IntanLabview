function [rate_map,rawRate_map,spatial_scale,SpksShuffle,oc_map,sinfo]=parpmap(Spks,Traj, msT,fstart,varargin)

p = inputParser;
p.addParamValue('animal','rat', @ischar);
p.addParamValue('shuffle', 0, @isnumeric);
p.addParamValue('shuffleType', 1, @isnumeric);
p.addParamValue('shufflen', 1000, @isnumeric);
p.addParamValue('verbose', 0, @isnumeric);
p.addParamValue('binside', 2.5, @isnumeric);

p.parse(varargin{:});
animal = p.Results.animal;
shuffle = p.Results.shuffle;
shuffleType = p.Results.shuffleType;
shuffleN=p.Results.shufflen;
verbose=p.Results.verbose;
binside=p.Results.binside;

year=2018;

shuffleLag=20000;%20sec
SpksShuffle=[];
BinWidthCm=2.5;%2.5cm
sinfo=[];

xl=3;
yl=4;

kHz=31.25;
ThS=2;%100cm/s
Linear=0;

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
end


switch lower(animal)
  case 'bird',
    cmPerPixel=120/Len;%120cm circle;
    Bin=BinWidthCm/cmPerPixel;
    spON=0;
  case 'sgbird',
    cmPerPixel=120/Len;%120cm circle;
    Bin=BinWidthCm/cmPerPixel;
    kHz=30;
    spON=0;          
  case 'birdsmall',
    cmPerPixel=80/Len;%120cm circle;
    Bin=BinWidthCm/cmPerPixel;
    spON=0;
  case 'birdsmalls',
    cmPerPixel=80/Len;%120cm circle;
    Bin=BinWidthCm/cmPerPixel;
    spON=1; 
    ThS=2.5;
  case 'birds',
    cmPerPixel=120/Len;%120cm circle;
    Bin=BinWidthCm/cmPerPixel;
    spON=1;
    ThS=2.5;
  case 'sgbirds',
    cmPerPixel=120/Len;%120cm circle;
    Bin=BinWidthCm/cmPerPixel;
    kHz=30;
    spON=1;
    ThS=2.5;
  case 'sg2birds',
    cmPerPixel=120/Len;%120cm circle;
    Bin=BinWidthCm/cmPerPixel;
    kHz=20;
    spON=1;
    ThS=2.5;
  case 'fishs',
    cmPerPixel=150/Len;%120cm circle;
    Bin=BinWidthCm/cmPerPixel;
    spON=1;
    ThS=2.5;%100cm/s
  case 'rat',
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


FPS=floor(1/(median(diff(msT))/1000));
fs_video=FPS;
msFPS=floor(1/FPS*1000);


[Spks,Traj,msT,speed,StartTraj,EndTraj]=SpkTraj(Spks,Traj,msT,fstart,spON,msFPS,kHz,cmPerPixel,ThS);

if Linear
    x=Traj;
    y=ones(size(Traj));
else
    x=ceil(Traj(:,xl));
    y=ceil(Traj(:,yl));
end


spatial_scale=cmPerPixel;
spk_x=[];
spk_y=[];

fieldShuffleFlag=0;


if shuffle
    
    if shuffleType==1 | shuffleType==2%SpikeTiming shuffle 
        beginSpks=fstart;%%
        endSpks=Spks(end);        
        
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
    elseif shuffleType==3%
        for i=1:shuffleN
            bootsamp=sort(randsample(length(Spks),length(Spks),true));
            SpksShuffle(i,:)=Spks(sort(bootsamp));
            BS(i,:)=bootsamp;
        end
    end
    
    if shuffleType==1 | shuffleType==3%spike timing shuffle or bootstrap
        for i=1:shuffleN

            Spks=SpksShuffle(i,:);
            spk_x=[];
            spk_y=[];
            for j=1:(length(msT)-1)
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

    elseif shuffleType==2%fieldShuffle
        rawRate_map=[];
        spk_x=[];
        spk_y=[];
        for j=1:(length(msT)-1)
            spkcnt=sum(Spks >= msT(j) & Spks < msT(j)+msFPS);
            if spkcnt
                for k=1:spkcnt
                    spk_x=[spk_x; Traj(j,xl)];
                    spk_y=[spk_y; Traj(j,yl)];
                end
            end
        end
        

        [rm,~,~,oc_map] = ratemap(spk_x,spk_y,x,y, ...
                                            spatial_scale,fs_video);            
        
        for i=1:shuffleN
            tmp=fieldShuffleEE(rm,oc_map,binside);
            [rate_map(i,:,:)] = tmp;
        end
        sinfo=[];
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
        plot(Traj(:,xl),Traj(:,yl),'k');
        hold on;
        scatter(spk_x,spk_y,40,'r','filled');
        axis equal off;
        set(gca,'xdir','normal');
        set(gca,'ydir','reverse');
        maxX=max(spk_x);
        maxY=max(spk_y);
        plot([maxX+10 maxX+10+1/cmPerPixel*20],[maxY+100 maxY+100],'k-');
    end
    

    [rm,~,~,oc_map] = ratemap(spk_x,spk_y,x,y, ...
                                            spatial_scale,fs_video);            

    rate_mapB=[];
    sinfo=[];
end
