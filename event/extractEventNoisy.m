%%
%%For nose-porking, events #1-3
%%[r,f]=extractEventNoisy(event(3,:),'threshold',2);
%%
%%for event #4 (treadmill)
%%[r,f]=extractEventNoisy(event(4,:),'Threshold',3,'pulsewidth',100);
%%
%%for event #5 (camera noiseless)
%%[r,f]=extractEventNoisy(event(5,:),'pulsewidth',100);
%%
function [r,f]=extractEventNoisy(event,varargin)
p = inputParser;
p.addParamValue('threshold', 3, @isnumeric);
p.addParamValue('pulsewidth', 10000, @isnumeric);

p.parse(varargin{:});
Th = p.Results.threshold;
pw = p.Results.pulsewidth;


se=smoothdata(event,'movmedian',100);
[r,f]=getTimes(se,Th);
id=find((f-r)>pw);
r=r(id);
f=f(id);
return;