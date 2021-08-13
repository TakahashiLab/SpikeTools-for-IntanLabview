function ensemble=addEnsemblePhase(ensemble,dlfp)
kHz=25;
binWidth=1;
sig_theta=BandpassFilter(dlfp,1000,[4 11]);

samplingRate=1000;
passband(1)=4;
passband(2)=11;
[wave,f,t,coh,wphases]=getWavelet(double(sig_theta),samplingRate,passband(1),passband(2),8,0);
[~,mIdx]=max(wave);%get index max power for each timepiont
pIdx=mIdx'+[0;size(f,2).*cumsum(ones(size(t,1)-1,1))];%converting to indices that will pick off single maxamp index from each of the freq-based phases at eacht timepoint
lfpphases=wphases(pIdx);%get phase of max amplitude wave at each timepoint
theta_phase = mod(lfpphases,2*pi);%covert to 0-2pi
                                  %rather than
theta_phase=rad2deg(theta_phase);
theta_phase=ceil(theta_phase./binWidth)';


for i=1:size(ensemble,1)
    ensemble{i,4}=theta_phase(floor(ensemble{i,3}/kHz));
end

return;