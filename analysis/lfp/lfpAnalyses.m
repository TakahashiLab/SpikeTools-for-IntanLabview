function lfpAnalyses(TrialT,dlfp,varargin)

p = inputParser;
p.addParamValue('dispmode', 'spect', @ischar);
p.addParamValue('channels', [1 1], @isvector);

p.parse(varargin{:});
dmode = p.Results.dispmode;
ch = p.Results.channels;


kHz=25;
states={'pre','peak','s','trough','s','rising','s','falling','post'};
maxFq=30;

srate_lfp=1000;

pre=[1 TrialT+1];
post=[TrialT size(dlfp,2)*kHz];

segments=ceil([pre' post']./kHz);

% Spectrogram parameters
WindowLength = 0.5; % Seconds
WindowLength = WindowLength*srate_lfp;
Overlap = round(0.45*srate_lfp);
NFFT = 2^13;


dlfpTarget=dlfp(ch(1),:);  
dlfpSource=dlfp(ch(2),:);  
%for i=1:size(segments,1)
c=1;
pSeq=[1 4 5 7 8 10 11 13 2];%plot
mSeq=[3 6 6 9 9 12 12 15 3];%mean
rSeq=[6 9 12 15];%reference plot
plvSeq=[1 3 4 5 6 7 8 9 2];

for i=1:size(segments,1)
%for i=1

  fprintf('seg #%d\n',i);

  seg=[segments(i,1):segments(i,2)];
  lfpTarget=dlfpTarget(seg);
  lfpSource=dlfpSource(seg);


  switch (lower(dmode))
    case 'spect',
      [SA, FA, TA, PA] = spectrogram(lfpTarget,...
                                     WindowLength,Overlap,NFFT, srate_lfp);

      FqPoint=max(find(FA<maxFq));
      subplot(5,3,pSeq(c));
      imagesc(TA,FA(1:FqPoint),PA(1:FqPoint,:));
      title(states{i});
      xlabel('time')
      ylabel('Frequency(Hz)');
      set(gca,'ydir','reverse');

      subplot(5,3,mSeq(c));
      hold on;
      plot(FA(1:FqPoint),mean(PA(1:FqPoint,:),2));

      [maxV,maxId]=max(mean(PA(1:FqPoint,:),2));
      maxId=FA(maxId);
      plot([maxId maxId],[0 maxV],'k');
      fprintf('%s:%f\n',states{i},maxV);
      for j=1:length(rSeq)
          if c==1 
              subplot(5,3,rSeq(j));
              hold on;
              plot(FA(1:FqPoint),mean(PA(1:FqPoint,:),2),'bo');
          elseif c==9
              subplot(5,3,rSeq(j));
              hold on;
              plot(FA(1:FqPoint),mean(PA(1:FqPoint,:),2),'go');
          end
      end
      c=c+1;
    case 'cohere',
      %  [Cxy_T,T,F]=cohereCore(lfpTarget,lfpSource,srate_lfp);

    case 'pac',
      maxFq=150;
      [xphaseT,xampT]=PAC(lfpTarget);
      [xphaseS,xampS]=PAC(lfpSource);
      out=getCFC(xphaseT(1:maxFq,:)',xampT');
      subplot(9,4,c)
      imagesc(out.MI);
      c=c+1;

      subplot(9,4,c)
      out=getCFC(xphaseS(1:maxFq,:)',xampS');
      imagesc(out.MI);
      c=c+1;

      subplot(9,4,c)
      out=getCFC(xphaseT(1:maxFq,:)',xampS');
      imagesc(out.MI);
      c=c+1;

      subplot(9,4,c)
      out=getCFC(xphaseS(1:maxFq,:)',xampT');
      imagesc(out.MI);
      c=c+1;

    case 'plvmod',
      subplot(5,2,plvSeq(c));
      plv_modindex_manager(lfpTarget',lfpSource',lfpTarget');
      title(states{i});
      c=c+1;
  end


end

colormap jet;



return;
%%%%%%%%%%%%%%%%%%%%%%
function [Cxy_T,T,F]=cohereCore(lfpTarget,lfpSource,srate_lfp)
%coherence
windowlength = 4*srate_lfp;
stepsize = round(srate_lfp/8);
Nwindow = (length(lfpTarget)-windowlength)/stepsize+1;
dt=1;
clear Cxy_T T
for nwin = 1:Nwindow
    win = (1:windowlength) + (nwin-1)* stepsize;
    [Cxy, F]=mscohere(lfpTarget(win),lfpSource(win),1*srate_lfp,...
                      [],10000,srate_lfp);
    Cxy_T(nwin,:)=Cxy;
    T(nwin) = median(win)*dt;
end

return;


