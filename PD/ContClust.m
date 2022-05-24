%edge: [start1 stop1;start2 stop2];
%edge(1,:)=>[start1 stop1]
%contEdge:continuous version of edge
function [edge,contEdge]=ContClust(seq)

[x,y]=size(seq);
if x>y
    seq=seq';
end

zeroOne=(diff(seq)==1);
edges=diff(zeroOne);

start=find(edges==1)+1;
stop=find(edges==-1)+1;

if zeroOne(1)==1
  start=[1 start];
end

if zeroOne(end)==1
  stop=[stop length(seq)];
end

edge=[start' stop'];

loop=size(edge,1);
contEdge=[];
for i=1:loop
  contEdge=[contEdge edge(i,1):edge(i,2)];
end

return;