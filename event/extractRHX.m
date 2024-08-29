function [NosePork,Treadmill,Pos,PosT]=extractRHX(event,t,Pos)

%nose-porking
for i=1:3
    [NosePork{i}]=extractEventNoisy(event(i,:),'threshold',2,'pulsewidth',100);
end

%%for event #4 (treadmill)
[Treadmill]=extractEventNoisy(event(4,:),'Threshold',3.2,'pulsewidth',1,'oversample',125,'swindow',10);

%%
%%for event #5 (camera noiseless)
[PosT]=extractEventNoisy(event(5,:),'pulsewidth',10);
delP=find(diff(t(PosT))*25>0.6);


Pos(delP+1,:)=[];


return;
