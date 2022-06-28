function [LD,LL,LR,P,D,pl]=PL(basename)

Suffix='.mat';
SuffixLen=size(Suffix,2)-1;

[path,name,ext]=fileparts(basename);
dataFolder=fullfile(path,name);
d=dir(fullfile(path,name));
loop=size(d,1);

LD=[];LL=[];LR=[];P=[];D=[];

possibleId=[];
for i=1:loop
    if length(d(i).name)>SuffixLen
        possibleId=[possibleId i];
    end
end
    
for i=possibleId
    if strncmp(d(i).name(end-SuffixLen:end),Suffix,SuffixLen)
        filename=fullfile(dataFolder,d(i).name);
        fprintf('loading %s\n',filename);
        load(filename,'LEDout','LFPoutL','LFPoutR','params','details');
        LD=[LD;LEDout];
        LL=[LL;LFPoutL];
        LR=[LR;LFPoutR];
        P=[P;params];
        D=[D;details];
    end
end

%beta
betaList=find(sum(str2mat(D{:,2})=='beta ',2)==5);
%theta
thetaList=find(sum(str2mat(D{:,2})=='theta',2)==5);
%gamma
gammaList=find(sum(str2mat(D{:,2})=='gamma',2)==5);

c={'r','g','b','k'};
fprintf('beta stimuli\n');

Data=LD;

pl=plotPL(betaList,Data,c{1});

fprintf('theta stimuli\n');
pl=plotPL(thetaList,Data,c{2});

fprintf('gamma stimuli\n');
pl=plotPL(gammaList,Data,c{3});
legend('beta','theta','gamma');

return;
%%%%%%%%%%%%%%
function pl=plotPL(List,Data,c)
bands{1}=1:4;%delta
bands{2}=5:12;%theta
bands{3}=13:30;%beta
bands{4}=31:61;%gamma
paraname={'Peak','P-S','Trough','T-S','Rising','R-S','Falling','F-S'};

pl=cell(4,2);
pl{1,2}='delta';
pl{2,2}='theta';
pl{3,2}='beta';
pl{4,2}='gamma';


data=Data(List,:,:);


for i=1:4
    for state=1:9
        pl{i,1}=[pl{i,1} mean(squeeze(data(:,state,bands{i})),2)];
    end
    [p,~,stats]=anova1(pl{i,1},[],'off');
    if p<0.05
        fprintf('slow phase: %s band \n',pl{i,2});
    end
    subplot(2,2,i);
    hold on;
    plot(mean(pl{i,1},1),c);
    title(['slow ' pl{i,2}]);

end



return;