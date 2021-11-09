function kkOuts=batchLoadKK(basename)

elecNums=1:8;

[path,name,ext]=fileparts(basename);

dataFolder=fullfile(path,name);
d=dir(fullfile(path,name));

for i=elecNums

    loadnameE=fullfile(dataFolder,[name 's' num2str(i) '.mat']);
    
    load(loadnameE,'kkOut');
    kkOuts{i}=kkOut;
end


return;
