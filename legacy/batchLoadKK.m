function kkOuts=batchLoadKK(basename)

elecNums=1:16;

[path,name,ext]=fileparts(basename);

dataFolder=fullfile(path,name);
d=dir(fullfile(path,name));

for i=elecNums

    loadnameE=fullfile(dataFolder,[name 's' num2str(i) '.mat']);
    if exist(loadnameE,'file')
        load(loadnameE,'kkOut');
        kkOuts{i}=kkOut;
    end
end


return;
