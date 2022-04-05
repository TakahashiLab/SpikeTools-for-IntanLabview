%%%
function [ensemble,Pos,PosT]=excludeDelay(ensemble,Up,Do,Pos,PosT,varargin)
p = inputParser;
p.addParamValue('binwidth', 100, @isnumeric);
p.addParamValue('samplingrate', 25, @isnumeric);
p.addParamValue('jitterpre', 0, @isnumeric);
p.addParamValue('jitterpost', 0, @isnumeric);

p.parse(varargin{:});
binWidth = p.Results.binwidth;
kHz=p.Results.samplingrate;
jitterPre=p.Results.jitterpre;
jitterPost=p.Results.jitterpost;


[x,y]=size(PosT);
if x<y
    msT=msT';
end

jitterPre=jitterPre*1000*kHz;
jitterPost=jitterPost*1000*kHz;

step=length(ensemble{1,1})/length(ensemble{1,3});

l=length(Up);
if l>length(Do);
    l=length(Do);
end

delIndT=[];
for j=1:size(ensemble,1);
    delInd=[];
    for i=l:-1:1
        orgSpks=ensemble{j,3};
        delInd=[delInd find(orgSpks>=Up(i) & orgSpks<=Do(i))];
        delIndT=[delIndT find(PosT>=Up(i) & PosT<=Do(i))'];
    end
    Seq=1:length(ensemble{j,3});
    Seq=Seq(delInd);
    ensemble{j,1}=getSpikesD(ensemble{j,1},step,Seq);
    ensemble{j,3}(delInd)=[];
end
for i=l:-1:1
    delIndT=[delIndT find(PosT>=Up(i) & PosT<=Do(i))'];
end

Pos(delIndT,:)=[];
PosT(delIndT)=[];

return;