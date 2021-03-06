function [speed,cmPerPixel]=createSpeed(Traj, msT,fstart,varargin)
verbose=0;
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


FPS=floor(1/(median(diff(msT))/1000));
fs_video=FPS;
msFPS=floor(1/FPS*1000);

if Linear
    dTraj=diff(Traj);
    ind=find(abs(dTraj)>max(Traj)/2);
    dTraj(ind)=dTraj(ind)-sign(dTraj(ind))*max(Traj);
    movement=sqrt(sum(dTraj.^2,2))*cmPerPixel;% 
else
    movement=sqrt(sum(diff(Traj(:,5:6)).^2,2))*cmPerPixel;% 
end

speed=movement./diff(msT);
speed=smooth(speed,FPS,'moving')*1000;
speed=[speed; speed(end)];

return;