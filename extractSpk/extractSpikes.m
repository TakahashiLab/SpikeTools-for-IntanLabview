function [tmp,tmps]=extractSpikes(filename,activeNums,step,Stimuli,filetype)

Th=50*(2^15-1)*5000*10^-6;%50uV 

if nargin<=3
  Stimuli=[];
  filetype='alvrc';
elseif nargin<=4
  filetype='ilvrc';
else
  
end

switch lower(filetype)
  case 'alvrc',
    Th=50*(2^15-1)*5000*10^-6;%50uV 
  case 'ilvrc',
    %    Th=300;%30uV 0.1uV
    Th=500;%50uV 0.1uV 500??? 200=20uV?
end

load(filename,'x');

%range=1500000; % 1 minutes at 25kHz

offset=0;
%Th=.5;
%Th=5.5;%%%%
%Th=7.5;

%estX=x(activeNums,1:range);
%estX=double(estX);
%stdV=-estimateNL(estX,Th);


stdV=-ones(12,1).*Th;
switch lower(filetype)
  case 'alvrc',
    x=-x(activeNums,:);%invert x for neuralynx Lynx-8 amplifiers
  case 'ilvrc',
    x=x(activeNums,:);
    stdV=-stdV;
end

stdV=stdV(activeNums);


[tmp,tmps]=extractSpDodeca4Bird(abs(x),stdV,step,x);
%[tmp,tmps]=extractSpDodeca(x,stdV,step);

if ~isempty(Stimuli)
    %  [tmp,tmps]=removeStimuli(tmp,tmps,step,Stimuli);
end


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kkOut]=normalSort(tmp,tmps,step)

NoE=size(tmp,1);
    
% feature extraction
fPCA=makePCA(tmp,size(tmp,1),step)';%'

% clustering with KlustaKwik
c=klustakwik(fPCA,1,10,1);

nC=max(c);
kkOut=[];
counter=1;
for i=1:nC
  Seq=find(c==i);
  out=getSpikesD(tmp,step,Seq);
  outs=tmps(Seq);
  if ~isempty(outs)
    kkOut{counter,1}=out;
    kkOut{counter,3}=outs;
    counter=counter+1;
  end
end

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%make principal component vectors 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fets] = makePCA(Spk, NoE, step)

nChannels = NoE; 
SpkSampls = step; 


nSpikes=size(Spk,2)/step;
Spk=reshape(Spk,[NoE,step,nSpikes]);

nPCA = 3; 
OnePCASet = 0; 


% go through channels
	for Channel=1:nChannels
		% calculate pca for this channel
		ChannelSpikes = squeeze(Spk(Channel, :, :))';%'
		CovMat = cov(ChannelSpikes);
        
        % calculate eigenvalues: why does it print by default?
        opts = struct('disp', 0);
		[ChannelPCs, d] = eigs(CovMat, eye(size(CovMat)),nPCA, 'LA', opts);
		
		% make it so PCs all have the same sign (hopefully...)
		ChannelPCs = ChannelPCs .* repmat(sign(ChannelPCs(13,:)), SpkSampls, 1);
		
		% now make feature vector
		ChannelFets = ChannelSpikes * ChannelPCs;
		
		% now add both of these to the output variables
		PCs(:,:,Channel) = ChannelPCs;
		Fets(:,(Channel-1)*nPCA + 1 : Channel*nPCA) = ChannelFets;
		
	end;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%add offsets
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tmps=addOffset(tmps,offset,step)

if offset~=0
  offset=offset-step;
end

tmps=tmps+offset;


return;


