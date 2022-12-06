%
function [DataE,ListE,Sub]=SPECTall(basename,varargin)
p = inputParser;
p.addParamValue('ref', 'LED', @ischar);
p.addParamValue('ledside', 'either', @ischar);
p.addParamValue('proc', 'individual', @ischar);
p.addParamValue('verbose', 1, @isnumeric);

p.parse(varargin{:});
ref = p.Results.ref;
ledside=p.Results.ledside;
proc=p.Results.proc;
verbose=p.Results.verbose;

global alpha;
alpha=0.05;

Suffix='.mat';
SuffixLen=size(Suffix,2)-1;

[path,name,ext]=fileparts(basename);
dataFolder=fullfile(path,name);
d=dir(fullfile(path,name));
loop=size(d,1);

LD=[];LL=[];LR=[];P=[];D=[];S=[];
AA=[];
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
    set(gca,'xtick',[1:16 ],'xticklabel',{'PD-stim','normal-stimu','PD-ref','normal-ref'});
    xlabel('\delta \theta \beta \gamma');

    set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
    ylabel('Feedback optogenetics');

    colormap bone;

  case 'individual',
    [resPost,resDuring]=processData(LD,LL,LR,D,S,ref,ledside,verbose);


    subplot(1,2,1);		
            image(resPost.*256);
            title(['Post stimulation (ref:' ref ',led:' ledside ')']);
            set(gca,'xtick',[1:4],'xticklabel',{'\delta','\theta','\beta','\gamma'});
            xlabel('power');
            set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
            ylabel('Feedback optogenetics');

            subplot(1,2,2);
            image(resDuring.*256);
            title(['During stimulation (ref:' ref ',led:' ledside ')']);
            set(gca,'xtick',[1:4],'xticklabel',{'\delta','\theta','\beta','\gamma'});
            xlabel('power');
            set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
            ylabel('Feedback optogenetics');
            colormap bone;



  case 'all',
    refs={'normal','pd','led'};
    ledsides={'normal','pd'};

        
    for i=1:length(refs)
	for j=1:length(ledsides)
            
            [resPost,resDuring,DataE{i,j},ListE{i,j},Sub{i,j}]=processData(LD,LL,LR,D,S,refs{i},ledsides{j},verbose);
            fprintf('ref:%s, ledside:%s\n',refs{i},ledsides{j});
            
            figure;
            subplot(1,2,1);		
            image(resPost.*256);
            title(['Post stimulation (ref:' refs{i} ',led:' ledsides{j} ')']);

            set(gca,'xtick',[1:4],'xticklabel',{'\delta','\theta','\beta','\gamma'});
            xlabel('power');
            set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
            ylabel('Feedback optogenetics');

            subplot(1,2,2);
            image(resDuring.*256);
            title(['During stimulation (ref:' refs{i} ',led:' ledsides{j} ')']);
            set(gca,'xtick',[1:4],'xticklabel',{'\delta','\theta','\beta','\gamma'});
            xlabel('power');
            set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
            ylabel('Feedback optogenetics');
            colormap bone;

        end
    end
end
return;
%%%%%%%%%%%%%%%

function [resPost,resDuring,Data,Lists,S]=processData(LD,LL,LR,D,S,ref,ledside,verbose)

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
S=S(Id,:);

%beta
betaList=find(sum(str2mat(D{:,2})=='beta ',2)==5);
%theta
thetaList=find(sum(str2mat(D{:,2})=='theta',2)==5);
%gamma
gammaList=find(sum(str2mat(D{:,2})=='gamma',2)==5);

Lists{1}=betaList;
Lists{2}=thetaList;
Lists{3}=gammaList;

resPost=[];
resDuring=[];

if verbose
    fprintf('###################################\n');
    fprintf('%%%%beta stimuli\n');
    fprintf('###################################\n');
end
[resP,resD]=plotPL(betaList,Data,c{1},S,verbose);
resPost=[resPost;resP];
resDuring=[resDuring;resD];

%save test.mat betaList Data
%return;
if verbose
    fprintf('###################################\n');
    fprintf('%%%%theta stimuli\n');
    fprintf('###################################\n');
end
[resP,resD]=plotPL(thetaList,Data,c{2},S,verbose);
resPost=[resPost;resP];
resDuring=[resDuring;resD];

if verbose
    fprintf('###################################\n');
    fprintf('%%%%gamma stimuli\n');
    fprintf('###################################\n');
end

[resP,resD]=plotPL(gammaList,Data,c{3},S,verbose);
resPost=[resPost;resP];
resDuring=[resDuring;resD];

return;
%%%%%%%%%%%%%%
function [resPost,resDuring]=plotPL(List,Data,c,S,verbose)

global alpha;

bands{1}=1:4;%delta
bands{2}=5:12;%theta
bands{3}=13:30;%beta
bands{4}=31:150;%gamma


pl=cell(4,2);
pl{1,2}='delta';
pl{2,2}='theta';
pl{3,2}='beta';
pl{4,2}='gamma';

state2{1,1}='left-left';
state2{1,2}='right-right';
state2{1,3}='left-right';
state2{1,4}='right-left';

paraname={'S','Peak','P-S','Trough','T-S','Rising','R-S','Falling','F-S'};

data=Data(List,:,:);
resPost=zeros(1,4);
resDuring=zeros(1,4);
Sub=S(List);

pacsize=size(data,3)/2;

xAxis=data(:,:,1:pacsize);
yAxis=data(:,:,1+pacsize:pacsize*2);

for i=1:4
    plc=[];
    sub=Sub;
    for state=1:9
        data=squeeze(yAxis(:,state,:));
        xa=squeeze(xAxis(1,state,:));

        id=find(xa>bands{i}(1) & xa<bands{i}(end));
        PLC=nanmean(data(:,id),2);
        plc=[plc PLC];
    end


    %delete outlier
    outlierId=[];
    for state=1:9
        if any(plc(:,state) > mean(plc(:,state))*10)
            outlierId=[outlierId; find(plc(:,state) > mean(plc(:,state))*10)];
        end
    end
    outlierId=unique(outlierId);
    plc(outlierId,:)=[];
    sub(outlierId,:)=[];

    
    if 0
    pn=paraname([1 3 5 7 9]);
    if plotStats(plc(:,[1 3 5 7 9]),pn,verbose)
        resPost(1,i)=1;
    end

    pn=paraname([1 2 4 6 8]);
    if plotStats(plc(:,[1 2 4 6 8]),pn,verbose)
        resDuring(1,i)=1;
    end
    else
        pn=paraname([1:9]);
        [resDuring(1,i),resPost(1,i)]=plotStats(plc(:,:),pn,sub,verbose);
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
function [resD,resP]=plotStats(plc,paran,sub,verbose)
global alpha;

resD=0;
resP=0;

if 0
    %[p,~,stats]=anova1(plc,[],'off');
    g1=reshape(repmat(1:size(plc,2),size(plc,1),1),1,prod(size(plc)));
    [p,~,stats]=anovan(reshape(plc,prod(size(plc)),1),{g1},'display','off');
    c=multcompare(stats,'CriticalValueType','dunnett','display','off');
    pv=c(:,6);
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
elseif 0
    t=table(plc(:,1),plc(:,2),plc(:,3),plc(:,4),plc(:,5));
    time=[1:5]';

    rm=fitrm(t,'Var1-Var5~1','withindesign',time);    
    [rtbl]=ranova(rm);
    p=rtbl.pValue;
    pv=multcompare(rm,'Time');
    c=[pv.Time_1(1:4) pv.Time_2(1:4)];
    pv=pv.pValue(1:4);
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

%if ~isempty(pp) & p < alpha
if ~isempty(pp) 
    for i=1:length(pp)
        if paran{c(pp(i),1)}=='S'
            if verbose
                fprintf('**********%s,p=%f\n',[paran{c(pp(i),1)} ' vs. ' paran{c(pp(i),2)}],pv(pp(i)));
            end
            if any(c(pp(i),2)==[2 4 6 8])
                resD=1;
            elseif any(c(pp(i),2)==[3 5 7 9])
                resP=1;
            end
        end
    end
end

return;

