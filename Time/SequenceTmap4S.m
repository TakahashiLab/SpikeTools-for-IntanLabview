function [seq,order]=SequenceTmap4S(ensemble,event,varargin)

p = inputParser;
p.addParamValue('l', 0, @isnumeric);
p.addParamValue('shuffle', 0, @isnumeric);
p.addParamValue('shufflen', 100, @isnumeric);
p.addParamValue('fstart', 0, @isnumeric);
p.addParamValue('jitterpre', 0, @isnumeric);
p.addParamValue('jitterpost', 0, @isnumeric);

p.parse(varargin{:});
L = p.Results.l;
shuffle=p.Results.shuffle;
shuffleN=p.Results.shufflen;
fstart=p.Results.fstart;
pre=p.Results.jitterpre;
post=p.Results.jitterpost;

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
        for j=1:shuffleN
            histR=plotTime(shuffleEns{i}(j,:),event,'verbose',0,'jitterpre',pre,'jitterpost',post,'smooth',1);
            org_seq(i,j,:)=histR;
        end
    else
        histR=plotTime(ensemble{i,3},event,'verbose',0,'jitterpre',pre,'jitterpost',post,'smooth',1);
        org_seq(i,:)=histR;
    end
end


if shuffle
    order=[];
    seq=org_seq1;
else
    [maxSeq,ind]=max(org_seq,[],2);
    [~,ind]=sort(ind);
    seq=org_seq./maxSeq;
    order=ind;
end

return;