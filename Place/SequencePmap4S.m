function [seq,order,oc_map]=SequencePmap4S(ensemble,LTraj,PosT,ThS,varargin)

p = inputParser;
p.addParamValue('outputmode', 'normalized', @ischar);
p.addParamValue('shuffle', 0, @isnumeric);
p.addParamValue('shufflen', 100, @isnumeric);
p.parse(varargin{:});
outputMode = p.Results.outputmode;
shuffle=p.Results.shuffle;
shuffleN=p.Results.shufflen;

if nargin<=3
    ThS=2.5;
end


loop=size(ensemble,1);
for i=1:loop
    fprintf('%d/%d\n',i,loop);

        if shuffle
            [seq(i,:,:)]=pmap(ensemble{i,3},LTraj,PosT,0,'animal','rat','speed',ThS,'verbose',0,'shuffle',shuffle,'shuffleN',shuffleN);
            oc_map=[];
        else
            [seq(i,:),~,~,oc_map]=pmap(ensemble{i,3},LTraj,PosT,0,'animal','rat','speed',ThS,'verbose',0);
        end

end

if shuffle
    outputMode='raw'
end
    
[maxSeq,ind]=max(seq,[],2);
[~,ind]=sort(ind);

switch (outputMode)
  case 'normalized',
    seq=seq./maxSeq;
  case 'raw',
end

order=ind;

return;