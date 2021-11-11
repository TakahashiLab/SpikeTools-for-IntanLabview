function [C,phi,f]=mtanalysiscCont(data1,data2,params)
Hz=25000;


%[C,phi,S12,S1,S2,f,zerosp,confC,phistd]=coherencyc(Data1,Data2,params);
[C,phi,S12,S1,S2,f]=coherencyc(data1,data2,params);
return;
for i=1:size(phistd,2)
  nanRange=find(phistd(:,i)>2);
  phi(nanRange,i)=NaN;
  C(nanRange,i)=NaN;
end


return;
