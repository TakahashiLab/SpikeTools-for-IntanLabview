function info=calcInfo(seq,oc_map)

if nargin==1
  PlaceMap=seq(1,:);
  oc_map=ones(size(PlaceMap));
end
for i=1:size(seq,1)
  PlaceMap=seq(i,:);
  [InfoPerSec,InfoPerSpk,InfoSparse]=calcInfoC(PlaceMap,oc_map);
  info(i)=InfoPerSpk;
end
return;
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [InfoPerSec,InfoPerSpk,InfoSparse]=calcInfoC(PlaceMap,OccupancyMap) 

tc=PlaceMap(:);
occ=OccupancyMap(:);
occ=occ/sum(occ);
f=sum(occ.*tc);
tc1=tc;
tc=tc/f;
ix=tc~=0;
SB=(occ(ix).*tc(ix)).*log2(tc(ix));
SB=sum(SB);

tc1=tc;
SS=(occ(ix).*tc1(ix)).*log2(tc(ix));
SS=sum(SS);

tau2=sum(tc1(ix).*occ(ix))^2;
SP=tau2/sum(occ(ix).*(tc1(ix).^2));
SP=sum(SP);

%%%%%%%%%%%%%%%%
taux=PlaceMap;
id=find(taux~=0 & ~isnan(taux));
taux=taux(id);
OccupancyMap=OccupancyMap(id);

px=OccupancyMap/sum(OccupancyMap);
%tau=sum(taux.*px);
tau=mean(taux);
%px

%meaning taui=PlaceMap;

if tau==0
  InfoPerSec=NaN;
  InfoPerSpk=NaN;
  InfoSparse=NaN;
  return;
end

%Info=(PlaceMap/tau);
tau2=sum(taux.*px)^2;

%tmp=taux.*log2(taux/tau).*px;
%surf(reshape(tmp,64,64)');

%InfoPerSec=sum(taux.*log2(taux/tau).*px);
%InfoPerSpk=sum((taux/tau).*log2(taux/tau).*px);
%InfoSparse=tau2/sum(px.*(taux.^2));
InfoPerSec=SS;
InfoPerSpk=SB;
InfoSparse=SP;

%taux
%(taux/tau)
%log2(taux/tau)

return;



