function sdat=filterAmp(x,beta,sampl,gpuFlag)
if nargin==1
  beta=0;
  sampl=20000;
  gpuFlag=1;
elseif nargin==3
    gpuFlag=1;
end

gpuOn=0;
if gpuDeviceCount & gpuFlag
    gpuOn=1;
end

switch(beta)
case 0%spikes
     bpFilt = designfilt('bandpassfir','FilterOrder',200, ...
         'CutoffFrequency1',800,'CutoffFrequency2',7500, ...
	 'SampleRate',sampl);
case 1%beta
     bpFilt = designfilt('bandpassfir','FilterOrder',200, ...
         'CutoffFrequency1',13,'CutoffFrequency2',30, ...
	 'SampleRate',sampl);
case 2%LFP
     bpFilt = designfilt('bandpassfir','FilterOrder',200, ...
         'CutoffFrequency1',0.1,'CutoffFrequency2',200, ...
	 'SampleRate',sampl);     
case 3 %beta2
     bpFilt = designfilt('bandpassiir','FilterOrder',2, ...
         'HalfPowerFrequency1',10,'HalfPowerFrequency2',30, ...
	 'SampleRate',sampl);     
case 4% theta
     bpFilt = designfilt('bandpassfir','FilterOrder',200, ...
         'CutoffFrequency1',4,'CutoffFrequency2',12, ...
	 'SampleRate',sampl);     
case 5% ripple
     bpFilt = designfilt('bandpassfir','FilterOrder',200, ...
         'CutoffFrequency1',80,'CutoffFrequency2',200, ...
	 'SampleRate',sampl);     
case 6% gamma
     bpFilt = designfilt('bandpassfir','FilterOrder',200, ...
         'CutoffFrequency1',30,'CutoffFrequency2',80, ...
                         'SampleRate',sampl);          
case 7%notch filter 50Hz
     bpFilt = designfilt('bandstopiir','FilterOrder',2, ...
         'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
                         'DesignMethod','butter','SampleRate',sampl);          
end
[b,a]=tf(bpFilt);

sdat=x;
for i=1:size(x,1)
    if gpuOn
        sdat(i,:)=gpuFiltFilt(b,a,x(i,:));
    else
        sdat(i,:)=filtfilt(bpFilt,double(x(i,:)));
    end
end

return;
     
