%
function [DataE,ListE,Sub]=Behvall(basename,varargin)
p = inputParser;
p.addParamValue('proc', 'individual', @ischar);
p.addParamValue('ledside', 'either', @ischar);
p.addParamValue('ref', 'LED', @ischar);
p.addParamValue('verbose', 1, @isnumeric);
p.addParamValue('alpha', 0.05, @isnumeric);

p.parse(varargin{:});

proc=p.Results.proc;
ref = p.Results.ref;
ledside=p.Results.ledside;
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
        load(filename,'LED','L','R','details');
        if ~isempty(LED)
            LD=[LD;LED];
            AA=[AA;LED];
        end
        LL=[LL;L];
        LR=[LR;R];        
        AA=[AA;L];
        AA=[AA;R];

        Lid=find(cell2mat(details(:,4))<=4);
        Rid=find(cell2mat(details(:,4))>=5);
        LEDid=find(cell2mat(details(:,3))==cell2mat(details(:,4)));
        
        if 1
            S=[S;ones(size(LED,1)+size(L,1)+size(R,1),1)*c];%subject
            D=[D;details(LEDid,:)];
            D=[D;details(setdiff(Lid,LEDid),:)];
            D=[D;details(setdiff(Rid,LEDid),:)];
        else
            switch(lower(ref))
              case 'normal',
                S=[S;ones(size(L,1),1)*c];%subject
                D=[D;details(setdiff(Lid,LEDid),:)];
              case 'pd',%
                S=[S;ones(size(R,1),1)*c];%subject
                D=[D;details(setdiff(Rid,LEDid),:)];
              case 'led',%
                S=[S;ones(size(LED,1),1)*c];%subject
                D=[D;details(LEDid,:)];
            end
        end
        c=c+1;
    end
end

switch lower(proc) 

    %The difference of optogenetics stimuli(led,left,right) is ignored
  case 'total',
    [resPost,resDuring]=processData2(AA,D,S,verbose);
    subplot(1,2,1);		
    image(resPost.*256);
    title(['Post stimulation']);
    set(gca,'xtick',[1:4 ],'xticklabel',{'PD-stim','normal-stimu','PD-ref','normal-ref'});
    xlabel('stim-ref');
    set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
    ylabel('Feedback optogenetics');


    subplot(1,2,2);
    image(resDuring.*256);
    title(['During stimulation']);
    set(gca,'xtick',[1:4 ],'xticklabel',{'PD-stim','normal-stimu','PD-ref','normal-ref'});
    xlabel('stim-ref');

    set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
    ylabel('Feedback optogenetics');

    colormap bone;

    %%% This option must specify refs and ledsides parameters 
  case 'individual',
    [resPost,resDuring]=processData(LD,LL,LR,AA,D,S,ref,ledside,verbose);

    subplot(1,2,1);		
            image(resPost.*256);
            title(['Post stimulation (ref:' ref ',led:' ledside ')']);
            %            set(gca,'xtick',[1],'xticklabel',{'behav'});
            set(gca,'xtick',1:4,'xticklabel',{'Peak','Trough','Rising','Falling'});
            xlabel('behavior');
            set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
            ylabel('Feedback optogenetics');

            subplot(1,2,2);
            image(resDuring.*256);
            title(['During stimulation (ref:' ref ',led:' ledside ')']);
            %            set(gca,'xtick',[1],'xticklabel',{'behav'});
            set(gca,'xtick',1:4,'xticklabel',{'Peak','Trough','Rising','Falling'});   

            xlabel('behavior');
            set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
            ylabel('Feedback optogenetics');
            colormap bone;


%%% This option processes and display all parameters of refs and ledsides at once. 
  case 'all',

    refs={'normal','pd','led'};
    ledsides={'normal','pd'};

        
    for i=1:length(refs)
	for j=1:length(ledsides)
            
            [resPost,resDuring,DataE{i,j},ListE{i,j},Sub{i,j}]=processData(LD,LL,LR,AA,D,S,refs{i},ledsides{j},verbose);
            fprintf('ref:%s, ledside:%s\n',refs{i},ledsides{j});
            
            figure;
            subplot(1,2,1);		
            image(resPost.*256);
            title(['Post stimulation (ref:' refs{i} ',led:' ledsides{j} ')']);


            set(gca,'xtick',1:4,'xticklabel',{'Peak','Trough','Rising','Falling'});
            xlabel('behavior');
            set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
            ylabel('Feedback optogenetics');

            subplot(1,2,2);
            image(resDuring.*256);
            title(['During stimulation (ref:' refs{i} ',led:' ledsides{j} ')']);

            set(gca,'xtick',1:4,'xticklabel',{'Peak','Trough','Rising','Falling'});
            xlabel('behavior');
            set(gca,'ytick',1:3,'yticklabel',{'\beta','\theta','\gamma'});
            ylabel('Feedback optogenetics');
            colormap bone;

        end
    end

end

return;
%%%%%%%%%%%%%%%
function [resPost,resDuring,Data,Lists,S]=processData(LD,LL,LR,AA,D,S,ref,ledside,verbose)

c={'r.','g.','b.','k.'};

Data=AA;
if 1
    switch(lower(ref))
      case 'normal',
        %        Data=LL;
        rid=find(cell2mat(D(:,4))<=4);
      case 'pd',%
                %        Data=LR;
        rid=find(cell2mat(D(:,4))>=5);
      case 'led',%
                 %        Data=LD;
        rid=find(cell2mat(D(:,3))==cell2mat(D(:,4)));
    end

    switch(lower(ledside))
      case 'either',
        Id=1:size(D,1);
      case 'normal',
        Id=find(cell2mat(D(:,3))<=4);
      case 'pd',
        Id=find(cell2mat(D(:,3))>=5);
    end

    Id=intersect(Id,rid);

    Data=Data(Id,:,:);
    D=D(Id,:);
    S=S(Id,:);

end

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
[resP,resD]=plotPL(betaList,Data,c{1},D,S,verbose);
resPost=[resPost;resP];
resDuring=[resDuring;resD];

%save test.mat betaList Data
%return;
if verbose
    fprintf('###################################\n');
    fprintf('%%%%theta stimuli\n');
    fprintf('###################################\n');
end
[resP,resD]=plotPL(thetaList,Data,c{2},D,S,verbose);
resPost=[resPost;resP];
resDuring=[resDuring;resD];

if verbose
    fprintf('###################################\n');
    fprintf('%%%%gamma stimuli\n');
    fprintf('###################################\n');
end

[resP,resD]=plotPL(gammaList,Data,c{3},D,S,verbose);
resPost=[resPost;resP];
resDuring=[resDuring;resD];

return;

%%%%%%%%%%%%%%%
function [resD,resP]=plotStats(plc,paran,det,sub,verbose)
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

    %sub=categorical(sub);
    LED=categorical(cell2mat(det(:,3))<=4);%true normal, false: pd
    REF=categorical(cell2mat(det(:,4))<=4);%true normal, false: pd
    t=table(plc(:,1),plc(:,2),plc(:,3),plc(:,4),plc(:,5),plc(:,6),plc(:,7),plc(:,8),plc(:,9),LED,REF);
    %    within=[1 2 2 2 2 2 2 2 2;1 2 3 2 3 2 3 2 3;1:9]';
    within=[1:9]';


    %rm=fitrm(t,'Var1-Var9~1','withindesign',within);    
    %rm=fitrm(t,'Var1-Var9~sub*LED*REF','withindesign',within);    
    rm=fitrm(t,'Var1-Var9~LED*REF','withindesign',within);    
    [rtbl]=ranova(rm,'withinmodel','Time');
    %p=rtbl.pValue(3);
    p=rtbl.pValue(1);
    pvLED=multcompare(rm,'Time','by','LED');
    pvREF=multcompare(rm,'Time','by','REF');

    %    c=[pv.Time_1(:) pv.Time_2(:)];
    %    pv=pv.pValue(:);
    %    pp=find(pv<alpha);
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

if 0
    if ~isempty(pp) 
        for i=1:length(pp)
            if paran{c(pp(i),1)}=='S'
                fprintf('**********%s,p=%f\n',[paran{c(pp(i),1)} ' vs. ' paran{c(pp(i),2)}],pv(pp(i)));
                if any(c(pp(i),2)==[2 4 6 8])
                    resD=1;
                    fprintf('%f\n',mean(plc(:,c(pp(i),[1 2]))));
                elseif any(c(pp(i),2)==[3 5 7 9])
                    resP=1;
                    fprintf('%f\n',mean(plc(:,c(pp(i),[1 2]))));
                end
            end
        end
    end
else

    resD=zeros(1,4);
    resP=zeros(1,4);
    %PD-stim
    pv=pvLED.pValue(pvLED.Time_1((pvLED.LED)=='false')==1);
    pp=find(pv<alpha);
    fprintf('####PD-stim\n');
    for i=1:length(pp)
        fprintf('**********S vs. %s,p=%f\n',paran{pp(i)+1},pv(pp(i)));
        if any(pp(i)+1==[2 4 6 8])
            resD(1,1)=1;
        elseif any(pp(i)+1==[3 5 7 9])
            resP(1,1)=1;
        end
    end

    %normal-stim
    pv=pvLED.pValue(pvLED.Time_1((pvLED.LED)=='true')==1);
    pp=find(pv<alpha);
    fprintf('####normal-stim\n');
    for i=1:length(pp)
        fprintf('**********S vs. %s,p=%f\n',paran{pp(i)+1},pv(pp(i)));
        if any(pp(i)+1==[2 4 6 8])
            resD(1,2)=1;
        elseif any(pp(i)+1==[3 5 7 9])
            resP(1,2)=1;
        end
    end

    %PD-REF
    pv=pvREF.pValue(pvREF.Time_1((pvREF.REF)=='false')==1);
    pp=find(pv<alpha);
    fprintf('####PD-REF\n');
    for i=1:length(pp)
        fprintf('**********S vs. %s,p=%f\n',paran{pp(i)+1},pv(pp(i)));
        if any(pp(i)+1==[2 4 6 8])
            resD(1,3)=1;
        elseif any(pp(i)+1==[3 5 7 9])
            resP(1,3)=1;
        end
    end

    %normal-REF
    pv=pvREF.pValue(pvREF.Time_1((pvREF.REF)=='true')==1);
    fprintf('####normal-REF\n');
    pp=find(pv<alpha);
    for i=1:length(pp)
        fprintf('**********S vs. %s,p=%f\n',paran{pp(i)+1},pv(pp(i)));
        if any(pp(i)+1==[2 4 6 8])
            resD(1,4)=1;
        elseif any(pp(i)+1==[3 5 7 9])
            resP(1,4)=1;
        end
    end

end

return;

