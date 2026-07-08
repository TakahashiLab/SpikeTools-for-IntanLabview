function batchPreProcessIntan(basename)

dirExist=0;
[path,name,ext]=fileparts(basename);

dataFolder=fullfile(path,name);
d=dir(fullfile(path,name));


loop=size(d,1);

for i=3:loop
  if d(i).isdir
    d2=dir(fullfile(dataFolder,d(i).name));

    loop2=size(d2,1);
    for j=3:loop2
      if d2(j).isdir
	d3=dir(fullfile(dataFolder,d2(j).name));
	fn=fullfile(dataFolder,d2(j).name);
	dirExist=1;
	fprintf('preProcessIlvrc(fn,3) layer 3\n%s\n\n',fn);
	preProcessIlvrc(fn,0);
	batchSelection(fn,1:8);
	makeEnsemble2(fn,1:8,8);
	makeEnsemble3(fn,1:8,8);
      end
    end
    
    if ~dirExist
      fn=fullfile(dataFolder,d(i).name);
      fprintf('preProcessIlvrc(fn,3) layer 2\n%s\n\n',fn);
      preProcessIlvrc(fn,0);
      batchSelection(fn,1:8);
      makeEnsemble2(fn,1:8,8);
      makeEnsemble3(fn,1:8,8);
      dirExist=0;
    end
  end    
end




return;