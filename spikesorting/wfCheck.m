function wfCheck(kkOuts,num)
close all;
for i=1:size(kkOuts{num},1)
    figure;
    dispDodeca(kkOuts{num},i,'limitlen',100);
end
return;