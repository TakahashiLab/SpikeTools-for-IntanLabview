function sdat=filterX(x,lowcut,highcut,sampl)
if nargin==3
  sampl=20000;
end

bpFilt = designfilt('bandpassfir','FilterOrder',200, ...
         'CutoffFrequency1',lowcut,'CutoffFrequency2',highcut, ...
	 'SampleRate',sampl);

sdat=filtfilt(bpFilt,x')';
     
return;
     
