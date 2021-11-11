function [Cs,phis,fs]=batchMTSc(Data,duration,tetPair)
Hz=25000;
params.Fs=1250;%2500
params.fpass=[0 130];
%params.tapers=[500 999];
%params.tapers=[60 119];
%params.tapers=[30 59];
params.tapers=[20 49];
params.err=[1 0.95];

step=Hz./params.Fs;

Data=Data(:,duration);
Data1=decimate(double(Data((tetPair(1)-1)*4+1,:)),step);
Data2=decimate(double(Data((tetPair(2)-1)*4+1,:)),step);


[Cs,phis,fs]=mtanalysiscCont(Data1,Data2,params);


return;