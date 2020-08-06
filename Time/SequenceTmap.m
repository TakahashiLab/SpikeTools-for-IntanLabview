function [seq1,seq2,order]=SequenceTmap(ensemble,event,nums)
if nargin==3
    order=1;%sequence order 
end
loop=size(ensemble,1);
wSize=10;

for i=1:loop
    fprintf('%d/%d\n',i,loop);
    histR=plotRaster(ensemble{i,3},event,nums,'verbose',0);
    seq1(i,:)=smooth(histR{1},wSize);
    seq2(i,:)=smooth(histR{2},wSize);
end

[maxSeq1,ind1]=max(seq1,[],2);
[~,ind1]=sort(ind1);
[maxSeq2,ind2]=max(seq2,[],2);
[~,ind2]=sort(ind2);

seq1=seq1./maxSeq1;
seq2=seq2./maxSeq2;
order=[ind1 ind2];
return;