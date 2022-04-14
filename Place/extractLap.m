%%[id1,id2]=extractLap(PosT,beNP{1},beNP{2});
function [id1,id2]=extractLap(PosT,event1,event2)

l=length(event1);
if length(event2) < l  
    l=length(event2);
end

id1=[];
for i=1:l
    id1=[id1; find(PosT>event1(i) & PosT<event2(i))];
end

id2=1:length(PosT);
id2=setdiff(id2',id1);

return;