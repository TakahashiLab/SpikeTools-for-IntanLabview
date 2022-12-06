%
function PACall(basename,varargin)
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



hemiS={'l-l','r-r','l-r','r-l'};
%phaseC={'\delta','\theta','\beta','gamma'};
phaseC={'\delta','\theta'};
%AmpC={'\delta','\theta','\beta','gamma'};
AmpC={'\beta','\gamma'};
lenPhase=length(phaseC);
lenAmp=length(AmpC);

X=sqrt(4*lenAmp*lenPhase);
Y=sqrt(4*lenAmp*lenPhase);
switch lower(proc) 
  case 'total',
    [resPost,resDuring]=processData2(AA,D,S,verbose);

    c=1;
    for i=1:4
        for j=1:lenAmp
            for k=1:lenPhase
                resPostD=resPost(:,1+(c-1)*4:c*4);
                subplot(X,Y,c);		
                c=c+1;
                image(resPostD.*256);
                title(['Post stimulation']);
                set(gca,'xtick',[1:4],'xticklabel',{'PD-stim','normal-stim','PD-ref','normal-ref'});
                xlabel([hemiS{i} phaseC{k} 'phase-' AmpC{j} 'amp']);

                set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
                ylabel('Feedback optogenetics');
            end
        end
    end
    colormap bone;

    figure;
    c=1;
    for i=1:4
        for j=1:lenAmp
            for k=1:lenPhase
                resDuringD=resDuring(:,1+(c-1)*4:c*4);
                subplot(X,Y,c);
                c=c+1;
                image(resDuringD.*256);
                title(['During stimulation']);
                set(gca,'xtick',[1:4],'xticklabel',{'PD-stim','normal-stim','PD-ref','normal-ref'});
                xlabel([hemiS{i} phaseC{k} 'phase-' AmpC{j} 'amp']);
                set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
                ylabel('Feedback optogenetics');
            end
        end
    end
    colormap bone;



  case 'individual',
    [resPost,resDuring]=processData(LD,LL,LR,D,S,ref,ledside,verbose);

    c=1;
    for k=1:4
        subplot(2,4,c);
        c=c+1;
        image(resPost([k 4+k 2*4+k],:).*256);
        title('Post stimulation');
        set(gca,'xtick',[1 5 9 13],'xticklabel',{'\delta','\theta','\beta','\gamma'});
        xlabel('phase-amp');
        set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
        ylabel('Feedback optogenetics');

        subplot(2,4,c);
        c=c+1;
        image(resDuring([k 4+k 2*4+k],:).*256);
        title('During stimulation');
        set(gca,'xtick',[1 5 9 13],'xticklabel',{'\delta','\theta','\beta','\gamma'});
        xlabel('phase-amp');
        set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
        ylabel('Feedback optogenetics');
    end
    colormap bone;

  case 'all',
    refs={'normal','pd','led'};
    ledsides={'normal','pd'};

    for i=1:length(refs)
	for j=1:length(ledsides)
            figure;
            [resPost,resDuring]=processData(LD,LL,LR,D,S,refs{i},ledsides{j},verbose);
            fprintf('ref:%s, ledside:%s\n',refs{i},ledsides{j});
	    c=1;

            for k=1:4
		subplot(2,4,c);		
                c=c+1;
                image(resPost([k 4+k 2*4+k],:).*256);
                %subplot(length(refs),length(ledsides)*2,c);
                %c=c+1;
                %image(resPost.*256);
                title(['Post stimulation (ref:' refs{i} ',led:' ledsides{j} ')']);
                set(gca,'xtick',[1 5 9 13],'xticklabel',{'\delta','\theta','\beta','\gamma'});
                xlabel('phase-amp');
                set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
                ylabel('Feedback optogenetics');

                %subplot(length(refs),length(ledsides)*2,c);
                subplot(2,4,c);
                c=c+1;
                image(resDuring([k 4+k 2*4+k],:).*256);
                %image(resDuring.*256);
                title(['During stimulation (ref:' refs{i} ',led:' ledsides{j} ')']);
                set(gca,'xtick',[1 5 9 13],'xticklabel',{'\delta','\theta','\beta','\gamma'});
                xlabel('phase-amp');
                set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
                ylabel('Feedback optogenetics');
            end
            colormap bone;
        end
    end

end

%legend('beta','theta','gamma');

return;
%%%%%%%%%%%%%%%
function [resPost,resDuring]=processData(LD,LL,LR,D,S,ref,ledside,verbose)

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
Sub=S(List);

resPost=zeros(4,16);
resDuring=zeros(4,16);

pacsize=150*150;
for k=1:4%%%state2l-l,r-r,l-r,r-l
    data2=data(:,:,1+pacsize*(k-1):pacsize*k);
    plc=[];
    for state=1:9
        data3=squeeze(data2(:,state,:));
        PLC=[];
        for j=1:4
            for i=1:4
                data3=reshape(data3,size(data2,1),150,150);
                PLC=[PLC mean(mean(data3(:,bands{i},bands{j}),3),2)];
            end
        end
        plc{state}=PLC;
    end

    if verbose
        fprintf('####pac %s\n',state2{1,k});
    end
    
    if 0
        for q=1:16
            PLC3=[];
            for z=[1 3 5 7 9]
                PLC3=[PLC3 plc{z}(:,q)];
            end

            if verbose
                fprintf('region=%d\n',q);
            end
            pn=paraname([1 3 5 7 9]);
            if plotStats(PLC3,pn,verbose)
                resPost(k,q)=1;
            end
        end

        if verbose
            fprintf('####pac %s\n',state2{1,k});
        end
        for q=1:16
            PLC3=[];
            for z=[1 2 4 6 8]
                PLC3=[PLC3 plc{z}(:,q)];
            end

            if verbose
                fprintf('region=%d\n',q);
            end

            pn=paraname([1 2 4 6 8]);
            if plotStats(PLC3,pn,verbose)
                resDuring(k,q)=1;
            end
        end
    else

        for q=1:16
            PLC3=[];
            for z=[1:9]
                PLC3=[PLC3 plc{z}(:,q)];
            end

            if verbose
                fprintf('region=%d\n',q);
            end

            pn=paraname([1:9]);
            [resDuring(k,q),resPost(k,q)]=plotStats(PLC3,pn,Sub,verbose);
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
function [resD,resP]=plotStats(plc,paran,sub,verbose)
global alpha;

resD=0;
resP=0;

if 0

    [p,~,stats]=anova1(plc,[],'off');
    c=multcompare(stats,'CriticalValueType','dunnett','display','off');
    pv=c(:,6);
    pp=find(pv<alpha);
else
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
end

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

