function seq=tetrodeMap(an,en)
seq=[];
for i=1:size(an,2)
    lenA=length(an{i});
    
    for j=1:size(en{i},1)
        lenA=lenA-(length(en{i}{j})-1);
    end
    seq=[seq ones(1,lenA)*i];
end

return