%
%LFPpower(LFPdata,PosT,assignments)
%assignments={opt-elec, opt-side, anti-opt-side}
%opt-elect: electrode which is optically stimulated
% opt-side:electrodes which is in the  hemisphere where the opt-electrode is inserted. 
% anti-opt-side: electrodes which is not in the hemisphere wheere the opt-electrode is inserted
% ex.assignments{1}=3; assingments{2}=[1:2 4]; assignments{3}=[5:8];
%
function mP=LFPpower(dlfp,PosT,ass)
if nargin<3  
  ass=[];
end
oFs=25000;
fs=1000;
mp=oFs/fs;

delta=[1 4];
theta=[4 12]; 
beta=[13 40];
gamma=[40 100];
movingwin=[0.5 0.016];
params.Fs=fs;
params.fpass=[0 100];
params.tapers=[5 9];

PosT=floor(PosT/mp)/1000;%frequency 25kHz -> 1kHz; & 1sec unit
loop2=size(PosT,2);


loop=size(dlfp,1);
Ss=[];
for i=1:loop
  [S,t,f]=mtspecgramc(dlfp(i,:),movingwin,params);
  Ss{i}=10*log10(S);
end


Ind=[];
for i=1:loop2
  [~,ind]=min(abs(t-PosT(i)));
  Ind=[Ind ind];
end


ranges{1}=delta;
ranges{2}=theta;
ranges{3}=beta;
ranges{4}=gamma;

for j=1:4
  mSs=[];
  for i=1:loop
    mS=mean(Ss{i}(:,find(f>ranges{j}(1) & f<ranges{j}(2))),2);
    mSs=[mSs mS];
  end


  bps=[];
  for tet=1:loop/4
    bps=[bps mean(mSs(:,1+(tet-1)*4:tet*4),2)];
  end  
  
  mP{j}=bps(Ind,:);
end


for j=1:4
  if isempty(ass)
    mP{j}=[mean(mP{j}(:,1:4),2) mean(mP{j}(:,5:8),2)];
  else
    mP{j}=[mP{j}(:,ass{1}) mean(mP{j}(:,ass{2}),2) mean(mP{j}(:,ass{3}),2)];
  end
end

return;