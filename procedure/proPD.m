%dn: 2018,2019,2020
function [LEDout,LFPoutL,LFPoutR,params,details]=proPD(dn,an,varargin)

p = inputParser;
p.addParamValue('method', 'plvmodm', @ischar);
p.addParamValue('region', 'all', @ischar);

p.parse(varargin{:});
method = p.Results.method;
region=p.Results.region;


    elecNum=4;%tetrode
    alpha=0.05;
    [homePath,dataPath]=PDdata(dn,an);
    TimingData='timing.mat';
    paraname={'Peak','P-S','Trough','T-S','Rising','R-S','Falling','F-S'};
    details=dataPath;

    loop=size(dataPath,1);
    params=cell(loop,2);

for i=1:loop
    %for i=1:2
    fprintf('analyzing %s\n',dataPath{i,1});
    [~,name]=fileparts(dataPath{i,1});
    loadname=fullfile(homePath,dataPath{i,1},[name 'LFP.mat']);
    load(loadname,'dlfp');
    loadname=fullfile(homePath,dataPath{i,1},TimingData);
    load(loadname,'segPara','TrialT');
    [~,mc,m]=plotRF(segPara,0);
    goodPara=find(mc(1:8,6)<alpha &  m(2:9)>m(1));
    badPara=find(mc(1:8,6)<alpha &  m(2:9)<m(1));
    params{i,1}=goodPara;
    params{i,2}=badPara;
    params{i,3}=dataPath{i,2};

    switch(lower(region))
      case  'all',
        candElec=1:8;
      case 'ctx',
        candElec=dataPath{i,5};
      case 'str',
        candElec=dataPath{i,6};
    end
    fprintf('LFP right\n');
    cand=setdiff(candElec,dataPath{i,3});
    cand1=min(cand(cand<=4));
    cand2=min(cand(cand>4));

    fprintf('LED stim\n');
    %%%LED simulus as a source channel
    if dataPath{i,3}<=4%Optical fiber was implanted in the left hemisphere
        if dataPath{i,4}<=4%LFP reference was on the left hemisphere
            LFPp=cand2;%We forced the LFP position to the right hemisphere
        else
            LFPp=dataPath{i,4};
        end
        ch=[1+(LFPp-1)*elecNum 1+(dataPath{i,3}-1)*elecNum];%LED as source
    else%Optical fiber was implanted in the right hemisphere
        if dataPath{i,4}<=4%LFP reference was on the left hemisphere
            LFPp=dataPath{i,4};
        else
            LFPp=cand1;%We forced the LFP position to the left hemisphere
        end
        ch=[1+(LFPp-1)*elecNum 1+(dataPath{i,3}-1)*elecNum];%LED as source
    end

    LEDout(i,:,:)=lfpAnalyses(TrialT,dlfp,'channels',ch,'dispmode',method,'verbose',0); 


    %%%%LED excluded, right=low frequency source
    fprintf('LFP right\n');
    ch=[1+(cand1-1)*elecNum 1+(cand2-1)*elecNum];%
    ch=sort(ch);

    LFPoutR(i,:,:)=lfpAnalyses(TrialT,dlfp,'channels',ch,'dispmode',method,'verbose',0); 

    fprintf('LFP left\n');
    %%%%LED excluded, left=low frequency source
    ch=[1+(cand1-1)*elecNum 1+(cand2-1)*elecNum];%LED as source
    ch=sort(ch,'descend');
    LFPoutL(i,:,:)=lfpAnalyses(TrialT,dlfp,'channels',ch,'dispmode',method,'verbose',0); 
    
end


return;