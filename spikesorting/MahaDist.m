function [IsolationDistance,LRatio]=MahaDist(Out,nums,ic)
if nargin<=2
  ic=4;
end
%target
%wf=double(Out{nums,1});%IC

wf=Out{nums,ic};%IC
wf=double(wf);
ts=Out{nums,3};

lenWf=length(wf);
lenTs=length(ts);

step=lenWf/lenTs;

targetFets=getFets(wf,step);
%size(targetFets)

%reference (not target)
len=size(Out,1);
wf=[];
for i=1:len
  if i~=nums
    wf=[wf Out{i,ic}];%IC
  end
end

wf=double(wf);


refFets=getFets(wf,step);

%D2=mahal([refFets; targetFets],targetFets);
%ascendD2=sort(D2(1:size(refFets,1)));

[xt,yt]=size(targetFets);
if xt <yt
IsolationDistance=NaN;
LRatio=NaN;
return;
end
D2=mahal(refFets,targetFets);
ascendD2=sort(D2);

if size(ascendD2,1) > size(targetFets,1)
  IsolationDistance=ascendD2(size(targetFets,1));
else
  IsolationDistance=NaN;
end

%LC=sum(1-chi2cdf(D2,8));
LC=sum(1-chi2cdf(D2,size(targetFets,2)));
LRatio=LC/size(targetFets,1);


return;
%%%%%%%%%%%%%%%%%%%%%%%%%  
function Fets=getFets(wf,step)
ewf=sqrt(wf.^2);
Energy=[];
nwf=wf;
lenTs=size(wf,2)/step;
for i=1:lenTs

  tmpEnergy=sum(ewf(:,1+(i-1)*step:i*step),2)/step;
  Energy=[Energy tmpEnergy];
  nwf(:,1+(i-1)*step:i*step)=nwf(:,1+(i-1)*step:i*step)./repmat(tmpEnergy,1,step);
end


PCA1=makePCAFirst(nwf,size(nwf,1),step);


Fets=[PCA1 Energy'];

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fets] = makePCAFirst(Spk, NoE, step)

nChannels = NoE; 
SpkSampls = step; 


nSpikes=size(Spk,2)/step;
Spk=reshape(Spk,[NoE,step,nSpikes]);

nPCA = 1; 
OnePCASet = 0; 


% go through channels
	for Channel=1:nChannels
		% calculate pca for this channel
		ChannelSpikes = squeeze(Spk(Channel, :, :))';%'
		CovMat = cov(ChannelSpikes);


        % calculate eigenvalues: why does it print by default?
	opts = struct('disp', 0);


		[ChannelPCs, d] = eigs(CovMat, eye(size(CovMat)),nPCA,'LA', opts);
		
		% make it so PCs all have the same sign (hopefully...)
		ChannelPCs = ChannelPCs .* repmat(sign(ChannelPCs(13,:)), SpkSampls, 1);
		
		% now make feature vector
		ChannelFets = ChannelSpikes * ChannelPCs;
		
		% now add both of these to the output variables
		PCs(:,:,Channel) = ChannelPCs;
		Fets(:,(Channel-1)*nPCA + 1 : Channel*nPCA) = ChannelFets;
		
	end;
return;