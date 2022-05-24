function [xphase,xamp]=PAC(unfilteredtrace)
sampl=1000;
step=25;

params.Fs=sampl;
params.fpass=[2 100];
params.trialave=1;
params.tapers=[5 9];
movingwin=[0.5 0.01];

%unfilteredtrace=decimate(double(unfilteredtrace),step);
lfp=unfilteredtrace;

startBand=1;
stopBand=60;%300
stepBand=2;
bCnt=1;
fprintf('%03d%% ', 0);
for band=startBand:stopBand
  fprintf('\b\b\b\b\b%03d%% ', floor(bCnt/(stopBand-startBand)*100));
  ft(bCnt,:)=filterX(lfp,band,band+stepBand,sampl);
  Hx(bCnt,:)=hilbert(ft(bCnt,:));
  xphase(bCnt,:)=angle(Hx(bCnt,:));
  xamp(bCnt,:)=abs(Hx(bCnt,:));
  bCnt=bCnt+1;
end
fprintf('\n');


return;
