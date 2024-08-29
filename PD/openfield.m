%EpsEdge, 1:ambulation, 2:immobility, 3:fine movement
%EpsDuration, 1:ambulation, 2:immobility, 3:fine movement
%rs: running speed
%Accel: Accelaration speed
function [EpsEdge,EpsDuration,rs,Accel]=openfield(Pos,PosT,aux,varargin)

p = inputParser;
p.addParamValue('hz', 25000, @isnumeric);
p.addParamValue('fps', 50, @isnumeric);
p.addParamValue('sth', 2, @isnumeric);
p.addParamValue('pixel2cm', 0.185, @isnumeric);
p.addParamValue('xy', [1 2], @isvector);

p.parse(varargin{:});
Hz = p.Results.hz;
fps=p.Results.fps;
sth=p.Results.sth;
pixel2cm = p.Results.pixel2cm;
xydim=p.Results.xy;

%Hz=25000;
%pixel2cm=0.185;%cm^2 = 1 pixel
%samplingRate=60;%50Hz; fps
%Th=2;%2cm/s Threshold rat?20
Th=sth;

ambulationPeriod=0.5;%0.5sec 
immobilityPeriod=1;%1sec
aux=double(aux);
accelTh=mean(aux(1,:))+std(aux(1,:));

if size(Pos,1)~=size(PosT,1)
  if size(Pos,1)<size(PosT,2)
    PosT(end)=[];
  end
  PosT=PosT';
end

traj=Pos(:,xydim).*pixel2cm;

trajD=diff(traj);
deltaTime=diff(PosT)./Hz;%sec

dist=sqrt(trajD(:,1).^2+trajD(:,2).^2);

runSpeed=(dist./deltaTime)';

rs=movingAverage(runSpeed,fps);

%ambulation
seqA=find(rs>Th);

[edges]=ContClust(seqA);		
ambulationEpisode=PosT(seqA(edges(find(diff(edges')>=ambulationPeriod*fps),:)));


if size(ambulationEpisode,2)==1 & size(ambulationEpisode,1)==2
  ambulationEpisode=ambulationEpisode';
end

edges4FineMov=edges(find(diff(edges')<ambulationPeriod*fps),:);

%immobility
seq=find(rs<=Th);

[edges]=ContClust(seq');

immobilityEpisode=PosT(seq(edges(find(diff(edges')>=immobilityPeriod*fps),:)));
loop=size(immobilityEpisode,1);
delEps=[];
for i=1:loop
  if max(aux(1,ceil(immobilityEpisode(i,1)/4):ceil(immobilityEpisode(i,2)/4)))>accelTh
%    delEps=[delEps i];
  end
end



fineMovEpisode=[PosT(seqA(edges4FineMov));PosT(seq(edges(find(diff(edges')<immobilityPeriod*fps),:)));immobilityEpisode(delEps,:)];


immobilityEpisode(delEps,:)=[];


%other fine movement

AmbEps=[];
loop=size(ambulationEpisode,1);
for i=1:loop
  AmbEps=[AmbEps ambulationEpisode(i,1):ambulationEpisode(i,2)];
end

ImmMovEps=[];
loop=size(immobilityEpisode,1);
for i=1:loop
  ImmMovEps=[ImmMovEps immobilityEpisode(i,1):immobilityEpisode(i,2)];
end
antiFineMov=[AmbEps ImmMovEps];


%FineMov=1:size(Pos,1);
FineMov=1:PosT(end);
FineMov=setdiff(FineMov,antiFineMov);

EpsEdge{1}=ambulationEpisode;
EpsEdge{2}=immobilityEpisode;
EpsEdge{3}=fineMovEpisode;
EpsDuration{1}=AmbEps;
EpsDuration{2}=ImmMovEps;
EpsDuration{3}=FineMov;

%%%%%%%
Comp=4;%aux sampling =1/4
Accel=cell(1,3);
for Eps=1:3
  loop=size(EpsEdge{Eps},1);
  for i=1:loop
    start=floor(EpsEdge{Eps}(i,1)/Comp);
    if start==0
      start=1;
    end
    stop=floor(EpsEdge{Eps}(i,2)/Comp);
    Accel{Eps}=[Accel{Eps} aux(:,start:stop)];
  end
end

return;