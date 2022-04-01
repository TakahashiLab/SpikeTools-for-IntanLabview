%beNP:the beginning edge of nose-porking
%tmr: treadmill speed rank
%tmSpeed: treadmill speed (pulse interval / samplingrate)
function [beNP,tmr,tmSpeed]=sawataniR1(NosePork,Treadmill)

beNP=cell(1,3);%beginning edge of nose-porking

%all porking points
allNp=cell2mat(NosePork);

%find the beginning of nose-porking
bp=min(allNp);

for i=1:3
    if any(NosePork{i}==bp)
        fpNum=i;%first porking number
        break;
    end
end


checkPoint=find(NosePork{fpNum}==bp);
checkTime=NosePork{fpNum}(checkPoint);
d=NosePork{fpNum}-checkTime;
NosePork{fpNum}(d<0)=[];
beNP{fpNum}=[beNP{fpNum} checkTime];

while 1
    %current porking number
    fpNum=mod(fpNum,3)+1;%because circular (1-3)
    d=(NosePork{fpNum}-checkTime);%
    NosePork{fpNum}(d<0)=[];    
    
    if length(d)<=1
        break;
    end
    [~,checkPoint]=min(d);
    
    checkTime=NosePork{fpNum}(checkPoint);
    beNP{fpNum}=[beNP{fpNum} checkTime];
end

%%%%%%%%%Treadmill speed rank
trialLen=length(beNP{2});
stateNum=3;
tmSpeed=[];
for trial=1:trialLen
    tmSpeed=[tmSpeed median(diff(Treadmill(Treadmill > beNP{2}(trial) & Treadmill < beNP{3}(trial))))];
end

factor=floor((max(tmSpeed)-min(tmSpeed))/stateNum);
speedState=ceil(tmSpeed/factor);
uSpeedState=unique(speedState);
for i=1:stateNum
    speedState(speedState==uSpeedState(i))=i;
end

tmr=speedState;%treadmill speed rank
return;
