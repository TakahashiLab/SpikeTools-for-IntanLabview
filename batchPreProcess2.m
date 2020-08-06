function batchPreProcess2(basename,MaxTetrode)

if nargin==1
    MaxTetrode=16;
end

%preprocess individual tetrodes one by one

for i=1:MaxTetrode
    fprintf('processing tetrode #%d\n',i);
    preProcessIlvrc2(basename,1,i);
end    
fprintf('processing event\n');
%process event signals
preProcessIlvrc2(basename,3);    
return;