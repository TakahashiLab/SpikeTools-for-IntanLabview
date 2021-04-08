% LMMの時には最後に1を追加
%Nsta=[11 21 31];,Nen=[20 30 40];とかにする

function [seq1,seq2,order,ori_seq1,ori_seq2]=SequenceTmap(ensemble,event,nums,jitter,PORK3,Nsta,Nen,varargin)
% p = inputParser;
% p.addParamValue('L', 0, @isnumeric);
% 
% L = p.Results.L;

p = inputParser;
p.addParamValue('l', 0, @isnumeric);
p.addParamValue('shuffle', 0, @isnumeric);
p.addParamValue('shufflen', 100, @isnumeric);
p.addParamValue('post', [], @ismatrix);
p.addParamValue('fstart', 0, @isnumeric);
p.parse(varargin{:});
L = p.Results.l;
shuffle=p.Results.shuffle;
shuffleN=p.Results.shufflen;
PosT=p.Results.post;
fstart=p.Results.fstart;


loop=size(ensemble,1);
wSize=10;

if shuffle
    if isempty(PosT)
        error('Please input PosT value\n');
        return;
    else
        for j=1:loop
            shuffleEns{j}=spikeShuffle(ensemble{j,3},PosT,fstart, ...
                                  'shuffleN',shuffleN);
        end
    end
end


for i=1:loop
    fprintf('%d/%d\n',i,loop);
    
    
    if shuffle
        if L == 0
            for j=1:shuffleN
                [~,histR]=plotRasterMM2(shuffleEns{j}(i,:),event,nums, ...
                                        PORK3,Nsta,Nen,'verbose', ...
                                        0,'jitter',jitter);

            end
        elseif L == 1
            for j=1:shuffleN
                [~,histR]=plotRasterLMM(shuffleEns{j}(i,:),event, ...
                                        nums,PORK3,Nsta,Nen,'verbose',0,'jitter',jitter);
            end
        end
        ori_seq1(i,j,:)=smooth(histR{1},wSize);
        ori_seq2(i,j,:)=smooth(histR{2},wSize);            
    else
        if L == 0
            [~,histR]=plotRasterMM2(ensemble{i,3},event,nums,PORK3,Nsta,Nen,'verbose',0,'jitter',jitter);
        elseif L == 1
            [~,histR]=plotRasterLMM(ensemble{i,3},event,nums,PORK3,Nsta,Nen,'verbose',0,'jitter',jitter);
        end
        ori_seq1(i,:)=smooth(histR{1},wSize);
        ori_seq2(i,:)=smooth(histR{2},wSize);
    end

end


if shuffle
    order=[];
    seq1=ori_seq1;
    seq2=ori_seq2;
else
    [maxSeq1,ind1]=max(ori_seq1,[],2);
    [~,ind1]=sort(ind1);
    seq1=ori_seq1./maxSeq1;
    [maxSeq2,ind2]=max(ori_seq2,[],2);
    [~,ind2]=sort(ind2);
    seq2=ori_seq2./maxSeq2;

    order=[ind1 ind2];
end

return;