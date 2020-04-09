function [tuning_curve, theta, wu2,ang_hd, mr,mvl,occupancy,num_spikes] = DirectionalTuningFnc(Spks,Traj, msT,fstart,varargin)
% Computes the head direction tuning curve. 
% 
% Arguments (see syntax (2)):
%
% binsize -> width of binsize in degrees
% Continuize -> merge overlapping epochs, and concatonate all data from
% root.epoch
%
% Returns
%
% tuning_curve -> NxMxO matrix where N is the number of bins as per
% 'binsize', O is the number of epochs and M is the number of cells
%
% theta -> vector of bin centers in degrees
% (1) [tuning_curve, theta] = root.DirectionalTuningFcn(cel);
% (2) [tuning_curve, theta] = root.DirectionalTuningFcn(cel, 'binsize', degrees, 'Continuize', 1);
% from https://github.com/hasselmonians/CMBHOME

binsize = .5;
theta=-180+binsize/2:binsize:180-binsize/2;
theta = theta(:);
    
verbose=0;
spON=0;
BinWidthCm=2.5;
kHz=31.25;
ThS=2;
yLen=max(Traj(:,2))-min(Traj(:,2));
xLen=max(Traj(:,1))-min(Traj(:,1));
Len=yLen;
if yLen < xLen
    Len=xLen;
end

if nargin==4
    cmPerPixel=120/Len;%120cm circle;
    Bin=BinWidthCm/cmPerPixel;
    spON=0;
    shuffle=0;
end

for i=1:2:(length(varargin)-1)
  % change the value of parameter
  switch lower(varargin{i})
    case 'animal',
      if strcmp(varargin{i+1},'bird')
          %Bin=12;
          cmPerPixel=120/Len;%120cm circle;
          Bin=BinWidthCm/cmPerPixel;
          spON=0;
      elseif strcmp(varargin{i+1},'sgbird')
          %Bin=12;
          kHz=30;
          cmPerPixel=120/Len;%120cm circle;
          Bin=BinWidthCm/cmPerPixel;
          spON=0;          
      elseif strcmp(varargin{i+1},'birdsmall')
          %Bin=12;
          cmPerPixel=80/Len;%120cm circle;
          Bin=BinWidthCm/cmPerPixel;
          spON=0;
      elseif strcmp(varargin{i+1},'birdsmalls')
          %Bin=12;
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
          %ThS=5;                     
      elseif strcmp(varargin{i+1},'rat')
          cmPerPixel=160/Len;%120cm circle;
          Bin=BinWidthCm/cmPerPixel;
          kHz=25;
          spON=1;
          ThS=2.5;          
      end
      shuffle=0;
    case 'shuffle'
      shuffle=1;
      SpksShuffle=varargin{i+1};
  end
end



FPS=floor(1/(median(diff(msT))/1000));
fs_video=FPS;
msFPS=floor(1/FPS*1000);

movement=sqrt(sum(diff(Traj(:,1:2)).^2,2))*cmPerPixel;% 
speed=movement./diff(msT);
speed=smooth(speed,FPS,'moving')*1000;

if spON
  Good= find(speed >= ThS );
  Traj=Traj(Good,:);
  msT=msT(Good);
end


Spks=ceil(Spks/kHz)+fstart;%msec
                           
Traj=Traj';
x=Traj(1,:);
y=Traj(2,:);
x2=Traj(3,:);
y2=Traj(4,:);
Traj=Traj';
headdir=atan2d(x-x2,y-y2)';

StartTraj=msT(1);
EndTraj=msT(end);

MRs=[];
WU2s=[];
MVLs=[];
if shuffle

    for i=1:size(SpksShuffle,1)

        Spks=SpksShuffle(i,:);
        Spks=Spks(find(Spks > StartTraj & Spks < EndTraj));
        spk_headdir=[];
        for j=1:(size(headdir,1)-1)
            spk_headdir=[spk_headdir ones(1,sum(Spks >= msT(j) & Spks < msT(j+1))).*headdir(j)];
        end

        [~,~,mr,wu2,mvl]=DirectionTuningCore(headdir,spk_headdir,binsize,fs_video,theta);
        MRs=[MRs mr];
        WU2s=[WU2s wu2];
        MVLs=[MVLs mvl];
    end
    tuning_curve=[];
    ang_hd=[];
    mr=MRs;
    wu2=WU2s;
    mvl=MVLs;
    occupancy=[];
    
else
    Spks=Spks(find(Spks > StartTraj & Spks < EndTraj));
    spk_headdir=[];
    for i=1:(size(headdir,1)-1)
        spk_headdir=[spk_headdir ones(1,sum(Spks >= msT(i) & Spks < msT(i+1))).*headdir(i)];
    end

    [tuning_curve,ang_hd,mr,wu2,mvl,occupancy,num_spikes]=DirectionTuningCore(headdir,spk_headdir,binsize,fs_video,theta);
end

return;
%%%%%%%%%%%%%%%%%%%%%%%%
function [tuning_curve,ang_hd,mr,wu2,mvl,angle_occupancy,num_spikes]=DirectionTuningCore(headdir,spk_headdir,binsize,fs_video,theta)

spk_headdir=spk_headdir';
%%%%%%%%%%%%%%%%%%%%%%%%head direction
angle_occupancy = DirectionalOccupancy(binsize, headdir,fs_video); 
num_spikes = hist(spk_headdir,theta);


window=29;

num_spikes=[fliplr(num_spikes(end-window+1:end)) num_spikes fliplr(num_spikes(1:window))];

angle_occupancy=[fliplr(angle_occupancy(end-window+1:end)) angle_occupancy fliplr(angle_occupancy(1:window))];
num_spikes = smooth(num_spikes,window,'moving')';
angle_occupancy = smooth(angle_occupancy,window,'moving')';

tuning_curve = num_spikes ./ angle_occupancy;

tuning_curve = tuning_curve(window+1:end-window);
angle_occupancy = angle_occupancy(window+1:end-window);
num_spikes = num_spikes(window+1:end-window);


[ang_hd, mr,mvl] = GetOtherStats(tuning_curve, theta); 
wu2 = WatsonsU2(spk_headdir, headdir);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ang_hd, mr,mvl] = GetOtherStats(tuning_curve, theta)

    tf = ~isnan(tuning_curve); % where we arent nan

    theta=theta*unitsratio('rad','deg');
    
    theta = theta(tf); % remove nans
    
    tuning_curve = tuning_curve(tf)';
    

    xs = tuning_curve.*cos(theta); % average 
    ys = tuning_curve.*sin(theta);
    

    ang_hd = atan2(mean(ys),mean(xs)); % mean direction

    mr = (cos(ang_hd)*sum(xs) + sin(ang_hd)*sum(ys)) / sum(tuning_curve); % mean resultant length
    
    mvl=sqrt(mean(ys)^2+mean(xs)^2);%mean vector length

return;