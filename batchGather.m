function batchGather(basename,elecNums,elecType)

if nargin==1
    elecNums=1:16;
    elecType='tetrode';
end

switch lower(elecType)
  case 'tetrode',
    numOfMicrowires=4;
  case 'stereotrode',
    numOfMicrowires=2;
end
    
%load all data
[path,name,ext]=fileparts(basename);
dataFolder=fullfile(path,name);

loop=length(elecNums);
loadnames=cell(1,loop);
xAll=[];

split=10;

for i=elecNums
    loadnames{i}=fullfile(dataFolder,[name 'e' num2str(i) '.mat']);
end

TM=[];
for j=1:split
    fprintf('Processing %d/%d\n',j,split);
    xAll=[];
    c=1;
    for i=elecNums
        load(loadnames{c},'x');
        c=c+1;
        if i==elecNums(1) 
            len=size(x,2);
            splitLen=floor(len/split);
        end
    
        if j==split
            xS=x(:,1+(j-1)*splitLen:end);
        else
            xS=x(:,1+(j-1)*splitLen:j*splitLen);
        end
        xAll=[xAll;xS];
    end
    tm = median(xAll,1);    
    TM=[TM tm];
end

c=1;
for i=elecNums
    load(loadnames{c},'x');
    x = bsxfun(@minus, x, TM); % subtract median of each time point    
    parsave(loadnames{c},'x',x);
    c=c+1;
end

return;