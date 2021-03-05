%%%%%%%%%%%
function speed=windowSpeed(Traj,msT,slideWin,cmPerPixel)

xlyl=1:2;%5:6

msFPS=median(diff(msT));
FPS=floor((1/msFPS)*1000);

if size(Traj,1)~=length(msT)
  seq=0:(length(Traj)-1);
  msT=seq*msFPS+msT(1);
  msT=msT';
end


speed=zeros(1,length(msT));
slideWin=5;%5
for k=1:slideWin
    if size(Traj,2)>1
        movement=sqrt(sum(diff(Traj(k:slideWin:end,xlyl)).^2,2))*cmPerPixel;
    else        
        %movement=sqrt(sum(diff(Traj(k:slideWin:end)).^2,2))*cmPerPixel;
        dTraj=diff(Traj(k:slideWin:end));
        ind=find(abs(dTraj)>max(Traj)/2);
        dTraj(ind)=dTraj(ind)-sign(dTraj(ind))*max(Traj);
        movement=sqrt(sum(dTraj.^2,2))*cmPerPixel;    
    end    
    wSpeed=movement./(msFPS*slideWin);
    rP=k:slideWin:length(speed);
    speed(rP(1:length(wSpeed)))=wSpeed;
end
speed=smooth(speed,FPS,'moving')*1000;

return;