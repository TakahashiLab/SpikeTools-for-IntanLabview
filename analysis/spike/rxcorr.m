function out=rxcorr(spike1,spike2,binWidth,dispRange,Hz,num1,num2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% an interface of the calculation of the normalized cross correlogram
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjust the order of 0.05msec range

dispRange=dispRange*Hz;

binWidth=binWidth*Hz; % adjust unit from msec to 0.05msec.
dispRange=floor(dispRange/binWidth);

spike1=floor(spike1/binWidth);
spike2=floor(spike2/binWidth);

[output]=rxccor(spike1,spike2,dispRange);

if (num1==num2)
   output(dispRange+1)=0;
end   

spike1=spike1(randperm(length(spike1)));

t=-dispRange:dispRange;
YPLIM=max(output)+1;
YMLIM=min(output);

titlename=makeTitle(num1,num2,binWidth/Hz);
xtick=[-dispRange,-dispRange*3/4,-dispRange/2,-dispRange/4,0,dispRange/4,dispRange/2,dispRange*3/4,dispRange];


out=cell(10,1);
out{1}=t;
out{2}=output;
out{3}=0;
out{4}=YPLIM;
out{5}=YMLIM;
out{6}=xtick;
out{7}=binWidth;
out{8}=0;
out{9}=titlename;
out{10}=dispRange;


clear spike;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Output=reshapeOutput(output,binWidth)
	Output=zeros(1,floor(length(output)/binWidth));	
	for i=1:floor(length(output)/binWidth);
		Output(i)=sum(output((i-1)*binWidth+1:i*binWidth));
	end	
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function name=makeTitle(num1,num2,bin)
name=sprintf('(%d vs. %d, bin=%2.1fmsec)',num1,num2,bin);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makeSpikeNumberName(spike1,spike2,count,dispRange,YPLIM);
	[basedn, spike1n, ext] = fileparts(spike1);
	[basedn, spike2n, ext] = fileparts(spike2);
	spike1n=strrep(spike1n,'_','-');
	spike2n=strrep(spike2n,'_','-');
	number=sprintf('%s=%d\n%s=%d',spike1n,count(1),spike2n,count(2));
	text(dispRange/2.2,YPLIM*4/5,number);
return;
