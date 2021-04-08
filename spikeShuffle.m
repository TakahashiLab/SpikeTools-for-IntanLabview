%
%rat only
%
function SpksShuffle=spikeShuffle(Spks,PosT,fstart,varargin)
p = inputParser;
p.addParamValue('khz', 25, @isnumeric);
p.addParamValue('shufflen', 100, @isnumeric);
p.parse(varargin{:});
kHz = p.Results.khz;
shuffleN=p.Results.shufflen;

shuffleLag=20000;%20sec


Spks=ceil(Spks/kHz)+fstart;%msec
PosT=PosT/kHz;

beginSpks=PosT(1);
endSpks=PosT(end);
entireLength=ceil(endSpks-beginSpks);

rp=randperm(entireLength-shuffleLag*2+1);

rp=rp(1:shuffleN);
orgRp=shuffleLag:(entireLength-shuffleLag);
rp=orgRp(rp);
lenOrgRp=length(orgRp);
lenSpks=length(Spks);
SpksShuffle=zeros(shuffleN,lenSpks);
for z=1:shuffleN
    spkss=Spks+rp(z);
    topSpks=find(spkss>=endSpks);
    spkss(topSpks)=spkss(topSpks)-endSpks+beginSpks;
    SpksShuffle(z,:)=sort(spkss);
end

return;
    