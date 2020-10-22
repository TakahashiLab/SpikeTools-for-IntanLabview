function [seq1,seq2,order]=SequenceTmap(ensemble,event,nums,jitter,splitNum)
if nargin==3
    order=1;%sequence order 
    jitter=0;
    splitNum=0;
elseif nargin==4
    splitNum=0;
end
loop=size(ensemble,1);
wSize=10;

for i=1:loop
    fprintf('%d/%d\n',i,loop);
    histR=plotRaster(ensemble{i,3},event,nums,'verbose',0,'jitter',jitter,'splitNum',splitNum);
    seq1(i,:)=smooth(histR{1},wSize);
%     
%     if size(histR,1)>1　三重野
      seq2(i,:)=smooth(histR{2},wSize);
%     end
end

[maxSeq1,ind1]=max(seq1,[],2);
[~,ind1]=sort(ind1);
seq1=seq1./maxSeq1;

% if size(histR,1)>1 9/3三重野
[maxSeq2,ind2]=max(seq2,[],2);
[~,ind2]=sort(ind2);
seq2=seq2./maxSeq2;
% else
%   seq2=[];
% ind2=[];
% end
order=[ind1 ind2];

return;