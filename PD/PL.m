%
function [LD,LL,LR,P,D]=PL(basename,varargin)
p = inputParser;
p.addParamValue('ref', 'LED', @ischar);
p.addParamValue('ledside', 'either', @ischar);
p.addParamValue('proc', 'individual', @ischar);
p.addParamValue('verbose', 1, @isnumeric);
p.addParamValue('alpha', 0.05, @isnumeric);

p.parse(varargin{:});
ref = p.Results.ref;
ledside=p.Results.ledside;
proc=p.Results.proc;
verbose=p.Results.verbose;

global alpha;
alpha=p.Results.alpha;


Suffix='.mat';
SuffixLen=size(Suffix,2)-1;

[path,name,ext]=fileparts(basename);
dataFolder=fullfile(path,name);
d=dir(fullfile(path,name));
loop=size(d,1);

LD=[];LL=[];LR=[];P=[];D=[];S=[];

possibleId=[];
for i=1:loop
    if length(d(i).name)>SuffixLen
        possibleId=[possibleId i];
    end
end
    
c=1;
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
        S=[S;ones(size(LEDout,1),1)*c];%subject
        c=c+1;
    end
end

switch(lower(ref))
  case 'normal',
    AA=LL;
  case 'pd',%
    AA=LR;
  case 'led',%
    AA=LD;
end


switch lower(proc) 
  case 'total',
    [resPost,resDuring]=processData2(AA,D,S,verbose);
    subplot(1,2,1);		
    image(resPost.*256);
    title(['Post stimulation']);
    set(gca,'xtick',[1:16],'xticklabel',{'PD-stim','normal-stimu','PD-ref','normal-ref'});
    xlabel('\delta \theta \beta \gamma');
    set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
    ylabel('Feedback optogenetics');


    subplot(1,2,2);
    image(resDuring.*256);
    title(['During stimulation']);
    set(gca,'xtick',[1:16],'xticklabel',{'PD-stim','normal-stimu','PD-ref','normal-ref'});
    xlabel('\delta \theta \beta \gamma');
    set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
    ylabel('Feedback optogenetics');

    colormap bone;

  case 'individual',
    [resPost,resDuring]=processData(LD,LL,LR,D,S,ref,ledside);
    subplot(1,2,1);
    image(resPost.*256);
    title('Post stimulation');
    set(gca,'xtick',1:4,'xticklabel',{'\delta','\theta','\beta','\gamma'});
    xlabel('slow oscillation');
    set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
    ylabel('Feedback optogenetics');

    subplot(1,2,2);
    image(resDuring.*256);
    title('During stimulation');
    set(gca,'xtick',1:4,'xticklabel',{'\delta','\theta','\beta','\gamma'});
    xlabel('slow oscillation');
    set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
    ylabel('Feedback optogenetics');

    colormap bone;
  case 'all',

    refs={'normal','pd','led'};
    ledsides={'normal','pd'};

    c=1;
    for i=1:length(refs)
        for j=1:length(ledsides)
            [resPost,resDuring]=processData(LD,LL,LR,D,S,refs{i},ledsides{j});
            subplot(length(refs),length(ledsides)*2,c);
            c=c+1;

            image(resPost.*256);
            title(['Post stimulation (ref:' refs{i} ',led:' ledsides{j} ')']);
            set(gca,'xtick',1:4,'xticklabel',{'\delta','\theta','\beta','\gamma'});
            xlabel('slow oscillation');
            set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
            ylabel('Feedback optogenetics');

            subplot(length(refs),length(ledsides)*2,c);
            c=c+1;
            image(resDuring.*256);
            title(['During stimulation (ref:' refs{i} ',led:' ledsides{j} ')']);
            set(gca,'xtick',1:4,'xticklabel',{'\delta','\theta','\beta','\gamma'});
            xlabel('slow oscillation');
            set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
            ylabel('Feedback optogenetics');
        end
    end
    colormap bone;
end

%legend('beta','theta','gamma');

return;
%%%%%%%%%%%%%%%
function [resPost,resDuring]=processData(LD,LL,LR,D,S,ref,ledside)

c={'r.','g.','b.','k.'};

switch(lower(ref))
  case 'normal',
    Data=LL;
  case 'pd',%
    Data=LR;
  case 'led',%
    Data=LD;
end

switch(lower(ledside))
  case 'either',
    Id=1:size(D,1);
  case 'normal',
    Id=find(cell2mat(D(:,3))<=4);
  case 'pd',
    Id=find(cell2mat(D(:,3))>=5);
end

Data=Data(Id,:,:);
D=D(Id,:);

%beta
betaList=find(sum(str2mat(D{:,2})=='beta ',2)==5);
%theta
thetaList=find(sum(str2mat(D{:,2})=='theta',2)==5);
%gamma
gammaList=find(sum(str2mat(D{:,2})=='gamma',2)==5);

resPost=[];
resDuring=[];
fprintf('###################################\n');
fprintf('%%%%beta stimuli\n');
fprintf('###################################\n');
[pl,resP,resD]=plotPL(betaList,Data,c{1},S);
resPost=[resPost;resP];
resDuring=[resDuring;resD];

fprintf('###################################\n');
fprintf('%%%%theta stimuli\n');
fprintf('###################################\n');
[pl,resP,resD]=plotPL(thetaList,Data,c{2},S);
resPost=[resPost;resP];
resDuring=[resDuring;resD];

fprintf('###################################\n');
fprintf('%%%%gamma stimuli\n');
fprintf('###################################\n');
[pl,resP,resD]=plotPL(gammaList,Data,c{3},S);
resPost=[resPost;resP];
resDuring=[resDuring;resD];

return;
%%%%%%%%%%%%%%
function [pl,resPost,resDuring]=plotPL(List,Data,c,S)

global alpha;

bands{1}=1:4;%delta
bands{2}=5:12;%theta
bands{3}=13:30;%beta
bands{4}=31:61;%gamma


pl=cell(4,2);
pl{1,2}='delta';
pl{2,2}='theta';
pl{3,2}='beta';
pl{4,2}='gamma';


paraname={'S','Peak','P-S','Trough','T-S','Rising','R-S','Falling','F-S'};

data=Data(List,:,:);
Sub=S(List);

resPost=zeros(1,4);
resDuring=zeros(1,4);
for i=1:4
    sub=Sub;
    for state=1:9
        pl{i,1}=[pl{i,1} mean(squeeze(data(:,state,bands{i})),2)];
    end

    plc=pl{i,1};    
    
    fprintf('####slow phase: %s band\n',pl{i,2});

    if 0
    PLC1=plc(:,1);
    PLC2=plc(:,[3 5 7 9]);
    PLC2=PLC2(:);

    p=plotStats2(PLC1,PLC2,pl{i,2},i,c,'post stimulus--');
    if p<alpha
        resPost(1,i)=1;
        PLC3=plc(:,[1 3 5 7 9]);
        pn=paraname([1 3 5 7 9]);
        plotStats(PLC3,pl{i,2},pn);
    end

    PLC1=plc(:,1);
    PLC2=plc(:,[2 4 6 8]);
    PLC2=PLC2(:);
    p=plotStats2(PLC1,PLC2,pl{i,2},i,c,'durig stimulus--');

    end

    %if p<alpha
    if 1
        resDuring(1,i)=1;
        if 0
            PLC3=plc(:,[1 2 4 6 8]);
            pn=paraname([1 2 4 6 8]);
            plotStats(PLC3,pl{i,2},pn);
        else
            PLC3=plc(:,[1:9]);
            pn=paraname([1:9]);
            [resDuring(1,i),resPost(1,i)]=plotStats(PLC3,pl{i,2},pn,sub);
        end
    end
end


return;
%%%%%%%%%%%%%%555
function p=plotStats2(plc1,plc2,plstr,i,c,stim)
global alpha;
    [p]=ranksum(plc1,plc2);
    if p<alpha
        fprintf('*******%s\n',stim);
    end
    if 0
    subplot(2,2,i);
    hold on;
    plot([median(plc1) median(plc2)],c);
    title(['slow ' plstr]);
    set(gca,'xtick',1:9,'xticklabel',paraname);
    end
return;
%%%%%%%%%%%%%%%
function [resD,resP]=plotStats(plc,plstr,paran,sub)
global alpha;
resD=0;
resP=0;

    if 0
        [p,~,stats]=anova1(plc,[],'off');
        if p<alpha
            c=multcompare(stats,'CriticalValueType','dunnett','display','off');
            pv=c(:,6);
            pp=find(pv<alpha);
        end
    elseif 0
        t=table(plc(:,1),plc(:,2),plc(:,3),plc(:,4),plc(:,5),plc(:,6),plc(:,7),plc(:,8),plc(:,9));
        %    within=[1 2 2 2 2 2 2 2 2;1 2 3 2 3 2 3 2 3;1:9]';
        within=[1:9]';
        rm=fitrm(t,'Var1-Var9~1','withindesign',within);    
        [rtbl]=ranova(rm,'withinmodel','Time');
        p=rtbl.pValue(3);
        pv=multcompare(rm,'Time');
        c=[pv.Time_1(:) pv.Time_2(:)];
        pv=pv.pValue(:);
        pp=find(pv<alpha);
    elseif 1
        t=table(plc(:,1),plc(:,2),plc(:,3),plc(:,4),plc(:,5),plc(:,6),plc(:,7),plc(:,8),plc(:,9),sub);
        %    within=[1 2 2 2 2 2 2 2 2;1 2 3 2 3 2 3 2 3;1:9]';
        within=[1:9]';
    
        %rm=fitrm(t,'Var1-Var9~1','withindesign',within);    
        rm=fitrm(t,'Var1-Var9~sub','withindesign',within);    
        [rtbl]=ranova(rm,'withinmodel','Time');
        %p=rtbl.pValue(3);
        p=rtbl.pValue(4);
        pv=multcompare(rm,'Time');
        c=[pv.Time_1(:) pv.Time_2(:)];
        pv=pv.pValue(:);
        pp=find(pv<alpha);
    else
            [p,tbl]=anova_rm({plc(:,2:5)' plc(:,1)'},'off');
            p=p(3);
            t=table(plc(:,1),plc(:,2),plc(:,3),plc(:,4),plc(:,5));
            time=[1:5]';
            rm=fitrm(t,'Var1-Var5~1','withindesign',time);    
            pv=multcompare(rm,'Time');
            c=[pv.Time_1(1:4) pv.Time_2(1:4)];
            pv=pv.pValue(1:4);
            pp=find(pv<alpha);
    end

    %        if ~isempty(pp) & p < alpha
    if ~isempty(pp) 
        for i=1:length(pp)
            if paran{c(pp(i),1)}=='S'
                fprintf('**********%s,p=%f\n',[paran{c(pp(i),1)} ' vs. ' paran{c(pp(i),2)}],pv(pp(i)));
                if any(c(pp(i),2)==[2 4 6 8])
                    resD=1;
                elseif any(c(pp(i),2)==[3 5 7 9])
                    resP=1;
                end
            end
        end
    end

    if 0
        subplot(2,2,i);
        hold on;
        plot(mean(plc,1),c);
        title(['slow ' plstr]);
        set(gca,'xtick',1:9,'xticklabel',paraname);
    end
 
return;

