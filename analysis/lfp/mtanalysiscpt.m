function [C,phi,f,confC,phistd,zerosp]=mtanalysiscpt(data,tetNum,unit,params,seq,step)
Hz=25000;
movingwin=[1 0.5];
unit=floor(unit./step);

loop=size(seq,2);
if loop==2
  loop=1;
end

duration=min(floor(diff(seq)/1000)*1000);
unitDuration=duration*step;

if 0 %trial average
  Data1=zeros(duration,loop);
  for i=1:loop
    Data1(:,i)=data(seq(i)+1:seq(i)+duration);
    %  spk=unit(unit>seq(i) & unit<seq(i)+unitDuration);
    spk=unit(unit>seq(i) & unit<seq(i)+duration);

    %  Spk(i).times=(double(spk')-seq(i))./Hz;
    Spk(i).times=(double(spk')-seq(i))./params.Fs;
  end
  [C,phi,S12,S1,S2,f,zerosp,confC,phistd]=coherencycpt(Data1,Spk,params,1);
else %segmentation every 1sec
  win=.5;%sec
  Data1=data(seq(1):seq(end)+duration)';
  spk=unit(unit>seq(1) & unit<seq(end)+duration);
  Spk.times=(double(spk')-seq(1))./params.Fs;%convert from kHz to Hz
  [C,phi,S12,S1,S2,f,zerosp,confC,phistd]=coherencysegcpt(Data1,Spk,win,params,1);
end


%[C,phi,S12,S1,S2,t,f,zerosp,confC,phistd]=cohgramcpt(Data1,Spk,movingwin,params,1);

for i=1:size(phistd,2)
  nanRange=find(phistd(:,i)>2);
%  phi(nanRange,i)=NaN;
%  C(nanRange,i)=NaN;
end


return;
