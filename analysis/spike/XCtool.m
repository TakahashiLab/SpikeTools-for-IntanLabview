%
%function XCTool(Out,varargin)
%
function XCtool(Out,varargin)
fontpara=0.5;
Hz=20;%20kHz sampling

%default
binWidth=0.1;
dispRange=2;
dataRange=1000;
Cal=1;
Corr=1;%correlation: raw(1) or normalize(2) 
%method='normalized';
method='raw';
period=[];
antiperiod=[];
printPara='off';
savePara='off';
disp='on';

for i=1:2:(length(varargin)-1)
    % change the value of parameter
    switch lower(varargin{i})
      case 'binwidth'
	binWidth=varargin{i+1};
	
      case 'disprange'
	dispRange=varargin{i+1};
	
      case 'hz'
	Hz=varargin{i+1};	
	
      case 'datarange'
	dataRange=varargin{i+1};	
	
      case 'combination'
	combi=varargin{i+1};
	if strcmp(combi,'all')
	  Cal=1;
	else  
	  Cal=2;
	  num1=combi(1);
	  num2=combi(2);
	end
	
      case 'ac'
	Cal=3;
	combi=varargin{i+1};
        num1=combi(1);
	num2=combi(2);
	
      case 'event' 
	Event=varargin{i+1};
	Corr=2;
	
      case 'method' 
	method=varargin{i+1};
	
      case 'period'
	period=varargin{i+1};
	Corr=1;
	
      case 'anti-period'
	antiperiod=varargin{i+1};
	Corr=1;	
	
      case 'print'
	printPara=varargin{i+1};
	
      case 'save'
	savePara=varargin{i+1};
	
      case 'display'
	disp=varargin{i+1};
	
      otherwise
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized parameter: ''' varargin{i} '''']);
    end;
end

     
len=size(Out,1);

switch Cal
  case 1
    counter=1;
    for i=1:len
      for j=i:len

	subplot(len,len,counter);
	spike1=double(sort(Out{i,3}));
	spike2=double(sort(Out{j,3}));


	switch Corr
	  case 1,
	    if ~isempty(period)
	      spike1=getSpPeriod(spike1,period);
	      spike2=getSpPeriod(spike2,period);
	    elseif ~isempty(antiperiod)
	      spike1=getSpAntiPeriod(spike1,period);
	      spike2=getSpAntiPeriod(spike2,period);
	    end

	    out=rxcorr(spike1,spike2,binWidth,dispRange,Hz,i,j);

	  case 2,
	    eventTime=double(Event);
	    out=nxcorr(spike1,spike2,eventTime,binWidth,dispRange,dataRange,method,'off',Hz,i,j);
	 end

	dispPattern(out,fontpara,Hz);
	counter=counter+1;
    
      end
      counter=counter+i;
    end	
    

  case 2%combination
    spike1=double(sort(Out{num1,3}));
    spike2=double(sort(Out{num2,3}));
    switch Corr
      case 1
	out=rxcorr(spike1,spike2,binWidth,dispRange,Hz,num1,num2);
	
	if strcmp(disp,'on')
	  scrz=get(0,'screensize');
	  figH=figure('Position',scrz);
	  set(figH,'PaperType','A4','PaperPositionMode','auto');
	  dispPattern(out,fontpara,Hz);	  
	end
	
	if ~strcmp(printPara,'off')
	  figurename=sprintf('%s-%d-%d.jpg',printPara,num1,num2);
	  eval(['print -djpeg99 ' figurename]);	  
	  close(figH);
	end
	
	if ~strcmp(savePara,'off')
	  savename=sprintf('%s-%d-%d.mat',savePara,num1,num2);
	  save(savename,'out');
	end
      case 2
	eventTime=double(Event);
	out=nxcorr(spike1,spike2,eventTime,binWidth,dispRange,dataRange,method,'off',Hz,num1,num2);

	if strcmp(disp,'on')
	  scrz=get(0,'screensize');
	  figH=figure('Position',scrz);
	  set(figH,'PaperType','A4','PaperPositionMode','auto');
	  dispPattern(out,fontpara,Hz);	  
	end
	
	if ~strcmp(printPara,'off')
	  figurename=sprintf('%s-%d-%d.jpg',printPara,num1,num2);
	  eval(['print -djpeg99 ' figurename]);	  
	  close(figH);
	end
	
	if ~strcmp(savePara,'off')
	  savename=sprintf('%s-%d-%d.mat',savePara,num1,num2);
	  save(savename,'out');
	end
	
    end
  case 3
    spike1=double(Out{num1,3});
    spike2=double(Out{num2,3});
    spike=sort([spike1 spike2]);
    out=axcorr(spike,binWidth,dispRange,Hz,num1);
    dispPattern(out,fontpara,Hz);
end


clear all;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spikeOut=getSpPeriod(spike,period)

len=size(period,1);
spikeOut=[];
for i=1:len
  pre=period(i,1);
  post=period(i,2);
  spikeOut=[spikeOut spike(find(spike > pre & spike < post))];
end

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spikeOut=getSpAntiPeriod(spike,period)

len=size(period,1);
spikeOut=[];
for i=1:len
  pre=period(i,1);
  post=period(i,2);
  spikeOut=[spikeOut spike(find(spike > pre & spike < post))];
end

spikeOut=setdiff(spike,spikeOut);

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispPattern(out,fontpara,Hz)
t=out{1};
outputNXC=out{2};
highConf=out{3};
YPLIM=out{4}+1;
YMLIM=out{5};
xtick=out{6};
binWidth=out{7};
count=out{8};
titlename=out{9};
dispRange=out{10};


if sum(outputNXC > highConf*1.5) >=1  | sum(outputNXC < -highConf*1.5) >=1
%  	plot(t,outputNXC,'k');
  bar(t,outputNXC,1);
else
%  	plot(t,outputNXC,'b');
  bar(t,outputNXC,1);
end


if ~isempty(highConf)
  hold on;
  plot(t,highConf);
  plot(t,-highConf);
end

title(titlename,'fontsize',20*fontpara);
axis([-dispRange,dispRange,YMLIM,YPLIM]);

xlabel('Time(msec)','fontsize',20*fontpara);
ylabel('spikes/bin','fontsize',20*fontpara);
set(gca,'Xtick',xtick,'XtickLabel',xtick*(binWidth/Hz),'fontsize',20*fontpara);

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeSpikeNumberName(spike1,spike2,count,dispRange,YPLIM,fontpara);
	[basedn, spike1n, ext] = fileparts(spike1);
	[basedn, spike2n, ext] = fileparts(spike2);
	spike1n=strrep(spike1n,'_','-');
	spike2n=strrep(spike2n,'_','-');
	number=sprintf('%s=%d\n%s=%d',spike1n,count(1),spike2n,count(2));
	text(dispRange/2.2,YPLIM*4/5,number,'fontsize',1*fontpara);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [numbers]=checkFromToParam(filename,method,fromto)

counter=[];
[basedn, basefn, ext] = fileparts(filename);
if isempty(method)
	checkfn=sprintf('%s_*.t',basefn);
	checkNode=dir(fullfile(basedn,checkfn));
	for i=1:length(checkNode)
		if isempty(findstr(checkNode(i).name,'a.t'))
			[dn,fn,ext]=fileparts(checkNode(i).name);
			[dummy,number]=strtok(fn,'_');
			number=strrep(number,'_','');
			counter=[counter str2num(number)];	
		end
	end	
else
	checkfn=sprintf('%s_*a.t',basefn);
	checkNode=dir(fullfile(basedn,checkfn));
	for i=1:length(checkNode)
		[dn,fn,ext]=fileparts(checkNode(i).name);
		[dummy,number]=strtok(fn,'_');
		number=strrep(number,'_','');
		number=strrep(number,'a','');
		counter=[counter str2num(number)];
	end	
end

if strcmp(lower(fromto),'all')==1
	numbers=sort(counter);		
elseif max(fromto) > max(counter)
	fprintf('The maximum number of your inputs(%d) is larger than that of neurons (%d)\n',max(fromto),max(counter));
	numbers=-1;
else
	numbers=fromto;
end
return;







