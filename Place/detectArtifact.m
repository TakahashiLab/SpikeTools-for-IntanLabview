%nw: index: normal spike waveform 
%da: index: detected artifact
function [nw,da]=detectArtifact(Spk,Traj,PosT,range)
binsize=2.5;

[~,spatial_scale]=pmap(Spk,Traj,PosT,0,'animal','rat', 'verbose', ...
                       0);

%Spk=Spk(Spk>=PosT(1) & Spk<=PosT(end));
Traj=Traj*spatial_scale./binsize;

indx=[];
for j=1:size(range,1)
    indx=[indx; find(Traj>range(j,1) & Traj<range(j,2))];
end


%plot(PosT,Traj,'.');
%hold on;

[edges]=ContClust2(indx',500);
%indx=indx(edges(3,1):edges(3,2));
%plot(PosT(indx),Traj(indx),'.r');

ED=PosT(indx(edges));

da=[];
for i=1:size(edges,1)
    da=[da find(Spk>ED(i,1) & Spk<ED(i,2))];

end

nw=setdiff(1:length(Spk),da);
return;