%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Spks,Traj,msT,speed,StartTraj,EndTraj]=SpkTraj(Spks,Traj,msT,fstart,spON,msFPS,kHz,cmPerPixel,ThS)

verbose=0;

if nargin==8
    ThS=2.5;
end

[x1,x2]=size(msT);
if x1<x2
    msT=msT';
end

if size(Traj,1)~=length(msT) 
  seq=0:(length(Traj)-1);
  msT=seq*msFPS+msT(1);
  msT=msT';
end

Spks=ceil(Spks/kHz)+fstart;%msec

slideWin=5;
speed=windowSpeed(Traj,msT,slideWin,cmPerPixel);

%avSpeed=calcAV(Traj,msT,slideWin,cmPerPixel);
%avSpeed=abs(avSpeed);

ThAV=40;%40

if verbose
    histogram(avSpeed,[0:200]);
end

%retention=calcRetention(Traj,msT,slideWin,cmPerPixel);
ThR=10;%2cm, 5
       %histogram(retention,[0:40]);

%histogram(retention);

%Th
ThS=2.5;

%analyzed if the speed is more than 5cm/s 
if spON
    Good= find(speed >= ThS );
    %    avGood= find(avSpeed <= ThAV);
    %    rGood=find(retention>=ThR);
    %    Good=intersect(Good,avGood);
    %%%Good=intersect(Good,rGood);
    Traj=Traj(Good,:);
    msT=msT(Good);
    speed=speed(Good);
end

StartTraj=msT(1);
EndTraj=msT(end);
Spks=Spks(find(Spks > StartTraj & Spks < EndTraj));

return;


