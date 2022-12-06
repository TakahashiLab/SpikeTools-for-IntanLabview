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
    SUB=categorical(sub);
    t=table(plc(:,1),plc(:,2),plc(:,3),plc(:,4),plc(:,5),plc(:,6),plc(:,7),plc(:,8),plc(:,9),LED,REF,SUB);
    %    within=[1 2 2 2 2 2 2 2 2;1 2 3 2 3 2 3 2 3;1:9]';
    within=[1:9;1 2 3 2 3 2 3 2 3;]';
    %rm=fitrm(t,'Var1-Var9~1','withindesign',within);    
    %rm=fitrm(t,'Var1-Var9~sub*LED*REF','withindesign',within);    
    rm=fitrm(t,'Var1-Var9~LED*REF+SUB','withindesign',within);    

    [rtbl]=ranova(rm,'withinmodel','w1+w2');
    fprintf('LED=%f,REF=%f\n',rtbl.pValue([2 3]));
    %p=rtbl.pValue(3);
    p=rtbl.pValue(1);
    pvLED=multcompare(rm,'w1','by','LED');
    pvREF=multcompare(rm,'w1','by','REF');
    pv2LED=multcompare(rm,'w2','by','LED');
    pv2REF=multcompare(rm,'w2','by','REF');


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

    stateC={'During','Post'};
    %PD-stim
    pv=pv2LED.pValue(pv2LED.w2_1==1 & (pv2LED.LED)=='false');
    pp=find(pv<alpha);
    if ~isempty(pp)
        fprintf('PD-Stim::');
        fprintf('%s=%f\n',stateC{pp},pv(pp));
    end

    pv=pv2LED.pValue(pv2LED.w2_1==1 & (pv2LED.LED)=='true');
    pp=find(pv<alpha);
    if ~isempty(pp)
        fprintf('normal-Stim::');
        fprintf('%s=%f\n',stateC{pp},pv(pp));
    end

    pv=pv2REF.pValue(pv2REF.w2_1==1 & (pv2REF.REF)=='false');
    pp=find(pv<alpha);
    if ~isempty(pp)
        fprintf('PD-REF::');
        fprintf('%s=%f\n',stateC{pp},pv(pp));
    end

    pv=pv2REF.pValue(pv2REF.w2_1==1 & (pv2REF.REF)=='true');
    pp=find(pv<alpha);
    if ~isempty(pp)
        fprintf('normal-REF::');
        fprintf('%s=%f\n',stateC{pp},pv(pp));
    end



    resD=zeros(1,4);
    resP=zeros(1,4);
    %PD-stim
    pv=pvLED.pValue(pvLED.w1_1==1 & (pvLED.LED)=='false');
    pp=find(pv<alpha);
    fprintf('####PD-stim\n');
    for i=1:length(pp)
        if any(pp(i)+1==[2 4 6 8])
            fprintf('Stim::');
            resD(1,1)=1;
        elseif any(pp(i)+1==[3 5 7 9])
            fprintf('Post::');
            resP(1,1)=1;
        end
        fprintf('**********S vs. %s,p=%f\n',paran{pp(i)+1},pv(pp(i)));
    end

    %normal-stim
    pv=pvLED.pValue(pvLED.w1_1==1 & (pvLED.LED)=='true');
    pp=find(pv<alpha);
    fprintf('####normal-stim\n');
    for i=1:length(pp)
        if any(pp(i)+1==[2 4 6 8])
            fprintf('Stim::');
            resD(1,2)=1;
        elseif any(pp(i)+1==[3 5 7 9])
            fprintf('Post::');
            resP(1,2)=1;
        end
        fprintf('**********S vs. %s,p=%f\n',paran{pp(i)+1},pv(pp(i)));
    end

    %PD-REF
    pv=pvREF.pValue(pvREF.w1_1==1 & (pvREF.REF)=='false');
    pp=find(pv<alpha);
    fprintf('####PD-REF\n');
    for i=1:length(pp)
        if any(pp(i)+1==[2 4 6 8])
            fprintf('Stim::');
            resD(1,3)=1;
        elseif any(pp(i)+1==[3 5 7 9])
            fprintf('Post::');
            resP(1,3)=1;
        end
        fprintf('**********S vs. %s,p=%f\n',paran{pp(i)+1},pv(pp(i)));
    end


    %normal-REF
    pv=pvREF.pValue(pvREF.w1_1==1 & (pvREF.REF)=='true');
    fprintf('####normal-REF\n');
    pp=find(pv<alpha);
    for i=1:length(pp)

        if any(pp(i)+1==[2 4 6 8])
            fprintf('Stim::');
            resD(1,4)=1;
        elseif any(pp(i)+1==[3 5 7 9])
            fprintf('Post::');
            resP(1,4)=1;
        end
        fprintf('**********S vs. %s,p=%f\n',paran{pp(i)+1},pv(pp(i)));
    end

end

return;

