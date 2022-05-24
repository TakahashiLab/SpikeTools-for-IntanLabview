function LFP=immobilityLFP(LFP,tetnum,PosT,immobilityEpisode)

region=[];
loop=size(immobilityEpisode,1);
for i=1:loop
  region=[region immobilityEpisode(i,1):immobilityEpisode(i,2)];
end

LFP=double(LFP((tetnum-1)*4+1,:));
LFP=LFP(region);

return;