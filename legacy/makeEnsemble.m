function ensemble=makeEnsemble(fn,an,en)
tetNum=16;

kkOuts=cell(1,tetNum);
for i=1:tetNum
    fprintf('load tet#%d\n',i);
    loadName=sprintf('%ss%d.mat',fn,i);
    load(loadName,'kkOut');
    kkOuts{i}=kkOut;
end
ensemble=makeEnsemble64(kkOuts,an,en);
return;