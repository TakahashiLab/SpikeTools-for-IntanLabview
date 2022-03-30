%
%
%filename: e?.mat rawdata
%          x?.mat extracted data  
%          s?.mat sorted data
%
%          mode: 0: load raw data->extraction->sorting
%                1: load raw data->extraction
%                2: load extracted data->sorting(ic)
%                3: load sorting(kk)->sorting(ic)
%                4: load extracted data->sorting(kk)
%                5: load raw data->extraction->sorting(kk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function batchSpikesIntan(basename,mode,elecNums,Th,suffix)

if nargin<3
    fprintf('input error\n');
    return;
elseif nargin==3
    Th=50;
    suffix='ilvrc';
elseif nargin==4
    suffix='ilvrc';
end


numOfMicrowires=4;%tetrode
numOfElectrodes=16;
step=32;
activeNums=numOfMicrowires;

if nargin<=1
  mode=0;
  elecNums=1:numOfElectrodes;
elseif nargin<=2
    elecNums=1:numOfElectrodes;
else
    numOfElectrodes=elecNums(end);
end


[path,name,ext]=fileparts(basename);

dataFolder=fullfile(path,name);
d=dir(fullfile(path,name));


%load activeNums
loadnameA=fullfile(dataFolder,['activeNums.mat']);
if exist(loadnameA,'file')
  load(loadnameA);
else
  for i=1:numOfElectrodes
    AllAN{i}=[1:numOfMicrowires];
  end
  fprintf('%s cannot be found\n',loadnameA);
end


Stimuli=[];%electrical stimulation artifact excluded in batchSpikes.m 

for i=elecNums

  activeNums=AllAN{i};
  if isempty(activeNums)
  else
    loadnameE=fullfile(dataFolder,[name 'e' num2str(i) '.mat']);
    
    %spike extraction
    if mode==0 | mode==1 | mode==5
      fprintf('Extracting spikes from electrode %d \n',i);
      savenameX=fullfile(dataFolder,[name 'x' num2str(i) '.mat']);

      if ~exist(savenameX,'file')

	[tmp,tmps]=extractSpikes(loadnameE,activeNums,step,Stimuli,suffix,Th);  
	save(savenameX,'tmp','tmps');
     else
	load(savenameX,'tmp','tmps');
      end
    elseif mode==2
      loadnameX=fullfile(dataFolder,[name 'x' num2str(i) '.mat']);
      load(loadnameX,'tmp','tmps');
    elseif mode==3
      loadnameX=fullfile(dataFolder,[name 'x' num2str(i) '.mat']);
      load(loadnameX,'tmp','tmps');
      loadnameS=fullfile(dataFolder,[name 's' num2str(i) '.mat']);
      load(loadnameS,'kkOut');
    elseif mode==4
      loadnameX=fullfile(dataFolder,[name 'x' num2str(i) '.mat']);
      load(loadnameX,'tmp','tmps');
    end
  
    if mode~=1
      fprintf('Sorting spikes from electrode %d \n',i);
      savenameS=fullfile(dataFolder,[name 's' num2str(i) '.mat']);
      
      if ~exist(savenameS,'file') | mode==3
	%spike sorting
	if mode==3
	  [OutputData,SampleOut,kkOut,OutputLWV]=ICSort(loadnameE,activeNums,tmp,tmps,'icasso','on','ck',kkOut,'isctest','on','filetype','ilvrc');
	elseif mode==4 | mode==5
	  [OutputData,SampleOut,kkOut]=ICSort(loadnameE,activeNums,tmp,tmps,'kkonly','on','filetype','ilvrc');
        else
	  [OutputData,SampleOut,kkOut,OutputLWV]=ICSort(loadnameE,activeNums,tmp,tmps,'icasso','on','isctest','on','filetype','ilvrc');
	  %      [OutputData,SampleOut,kkOut,OutputLWV]=ICSort(loadnameE,activeNums,tmp,tmps,'icasso','off','isctest','on');
	end
	save(savenameS,'OutputData','SampleOut','kkOut','-v7.3');
      end
    end
  end
end

%delete(pp);
return;