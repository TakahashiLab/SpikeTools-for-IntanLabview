% seq0:元のseqもしくはどちらかのseq,seq1:前半,seq2:後半
%前半と後半の相関が0.1%以下なら

function [placeCell,order,seq,ind1,ind2,seq1,seq2]=seqCorr(seqP1,seqP2,seqP0,seqSP1,seqSP2,varargin)
% function [placeCell]=seqCorr(seqP1,seqP2)
alpha=0.01;
debug=0;

percentile=.95;
l1=size(seqP1,2);
l2=size(seqP2,2);
l3=size(seqP0,2);
l=min([l1 l2 l3]);
l=min([l1 l2]);
seqP1=seqP1(:,1:l);
seqP2=seqP2(:,1:l);
seqP0=seqP0(:,1:l);

cellNum=size(seqP1,1);

placeCell=[];
for i=1:cellNum
    ind1=find(~isnan(seqP1(i,:)));
    ind2=find(~isnan(seqP2(i,:)));
    ind=intersect(ind1,ind2);
    if ~isempty(ind)
        [CC,pval]=corr(seqP1(i,ind)',seqP2(i,ind)','Type','Pearson');
    else
        CC=NaN;
        pval=NaN;
    end
    if debug
        plot(seqP1(i,ind));
        hold on;
        plot(seqP2(i,ind));
        CC
        pval
    end
    
    sCC=[];
    for j=1:size(seqSP1,2)
        ind1=find(~isnan(seqSP1(i,j,:)));
        ind2=find(~isnan(seqSP2(i,j,:)));
        ind=intersect(ind1,ind2);
        scc=corr(squeeze(seqSP1(i,j,ind)),squeeze(seqSP2(i,j,ind)),'Type','Pearson');
        sCC=[sCC scc];  
    end
    sCC=sort(sCC);
    cutoff=sCC(length(sCC)*percentile);
    %     if pval<alpha 
    %if CC>0.3
    if CC>cutoff
        placeCell=[placeCell i];
    end
end

if nargin == 3
    seq=seqP0(placeCell,:);
    [maxSeq,ind]=max(seq,[],2);
    [~,ind]=sort(ind);
    seq=seq./maxSeq;
    order=ind;
    ind1=[];
    ind2=[];
    seq1=[];
    seq2=[];
elseif nargin == 4
    seq1=seqP1(placeCell,:);
    [maxSeq1,ind1]=max(seq1,[],2);
    [~,index1]=sort(ind1);
    seq1=seq1./maxSeq1;
    order1=index1;
    order{1}=order1;
    
    
    seq2=seqP2(placeCell,:);
    [maxSeq2,ind2]=max(seq2,[],2);
    [~,index2]=sort(ind2);
    seq2=seq2./maxSeq2;
    order2=index2;
    order{2}=order2;
    
    seq=[];
end
    

return;