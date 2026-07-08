function dlfp=downlfp(lfp,downN,kHz)

if nargin==2
  kHz=25;
end

%dlfp=zeros(size(lfp,1),ceil(size(lfp,2)/downN));
movingwin=[0.3 0.1];
params.Fs=kHz*1000;
params.tapers=[10 19];
params.fpass=[0 100];

for i=1:size(lfp,1)
      fprintf('%d/%d channel\n',i,size(lfp,1));
%  rlfp=rmlinesmovingwinc(double(lfp(i,:)),movingwin,10,params,[],[],60);
%  dlfp(i,:)=decimate(rlfp,downN);
  dlfp(i,:)=decimate(double(lfp(i,:)),downN);
end

[dlfp]=locdetrend(dlfp',1000,movingwin)';

return;