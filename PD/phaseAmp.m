function lfps=phaseAmp(dlfp,band,EpsEdge,Eps)
if nargin==2
  EpsEdge=[];
  Eps=[];
end

sampl=1000;
window=1000;%1sec
params.Fs=sampl;
params.fpass=[5 100];
params.trialave=1;
params.taper=[5 9];
%movingwin=[0.5 0.01];

bands=[
  3 6;%delta
  4 12;%theta
  13 30;%beta
  30 60;%lowgamma
  60 120];%highgamma

if ~isempty(EpsEdge)
  cEps=floor(EpsEdge{Eps}./25);%%%
  duration=[];
  for j=1:size(cEps,1)
    duration=[duration cEps(j,1):cEps(j,2)];
  end

  dlfp=dlfp(duration);
end

  fdlfp=filterX(dlfp,bands(band,1),bands(band,2),sampl);
  Hx=hilbert(fdlfp);
  xp=angle(Hx);
  xa=abs(Hx);
  
  
  [~,zP]=findpeaks(abs(xp),'minpeakWidth',125,'minpeakheight',3);

  plfp=sqrt(dlfp.^2);
  th=mean(plfp)+3*std(plfp);
%  window=ceil(1/bands(band,2)*sampl)*1;
  lfps=[];
  for i=1:length(zP)
    
    if (zP(i)-window)<0 | (zP(i)+window)>size(dlfp,2)
    else
      checkD=plfp(zP(i)-window:zP(i)+window);
%      if sum(checkD>th) >0
	lfps=[lfps;dlfp(zP(i)-window:zP(i)+window)];
%      end
    end 
  end
  
   lfps=lfps';

%movingwin=[window/50/sampl window/50/sampl/10]; 
movingwin=[0.5 0.05];
  [S,t,f]=mtspecgramc(lfps,movingwin,params);
  imagesc(t,f,10*log10(S)');
  set(gca,'ydir','normal');
  colormap(jet);


return;
