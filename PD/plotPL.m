
%%%%%%%%%%%%%%
function [resPost,resDuring]=plotPL(List,Data,c,D,S,verbose)

global alpha;

bands{1}=1:4;%delta
bands{2}=5:12;%theta
bands{3}=13:30;%beta
bands{4}=31:60;%gamma 150 for SPECT?


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
%phaseC={'\delta','\theta','\beta','gamma'};
phaseC={'\delta','\theta'};
%AmpC={'\delta','\theta','\beta','gamma'};
AmpC={'\beta','gamma'};
lenPhase=length(phaseC);
lenAmp=length(AmpC);

sub=S(List);
det=D(List,:);
pn=paraname([1:9]);
if length(size(Data))==3
    data=Data(List,:,:);
    resDuring=[];
    resPost=[];

    if size(data,3)>100*150
        pacsize=150*150;
        for k=1:4%%%state2l-l,r-r,l-r,r-l
            data2=data(:,:,1+pacsize*(k-1):pacsize*k);
            plc=[];
            for state=1:9
                data3=squeeze(data2(:,state,:));
                PLC=[];
                for j=1:lenAmp
                    for i=1:lenPhase
                        data3=reshape(data3,size(data2,1),150,150);
                        PLC=[PLC mean(mean(data3(:,bands{i},bands{j}),3),2)];
                    end
                end
                plc{state}=PLC;
            end


            q=1;
            for j=1:lenAmp
                for i=1:lenPhase
                    PLC3=[];
                    for z=[1:9]
                        PLC3=[PLC3 plc{z}(:,q)];
                    end
                    q=q+1;
                    fprintf('%s phase- %s Amp\n',phaseC{i},AmpC{j});
                    pn=paraname([1:9]);
                    [resD,resP]=plotStats(PLC3,pn,det,sub,verbose);
                    resDuring=[resDuring resD];
                    resPost=[resPost resP];
                end
            end
        end
    elseif size(data,3)>2000
        pacsize=size(data,3)/2;
        xAxis=data(:,:,1:pacsize);
        yAxis=data(:,:,1+pacsize:pacsize*2);

        for i=1:4
            plc=[];
            for state=1:9
                data=squeeze(yAxis(:,state,:));
                xa=squeeze(xAxis(1,state,:));

                id=find(xa>bands{i}(1) & xa<bands{i}(end));
                PLC=nanmean(data(:,id),2);
                plc=[plc PLC];
            end
            [resD,resP]=plotStats(plc,pn,det,sub,verbose);
            resDuring=[resDuring resD];
            resPost=[resPost resP];
        end
    else
        
        for i=1:4
            plc=[];
            for state=1:9
                plc=[plc mean(squeeze(data(:,state,bands{i})),2)];
            end
            fprintf('########################%s####################\n',pl{i,2});
            [resD,resP]=plotStats(plc,pn,det,sub,verbose);
            resDuring=[resDuring resD];
            resPost=[resPost resP];
        end
    end
else
    plc=Data(List,:);
    [resDuring,resPost]=plotStats(plc,pn,det,sub,verbose);
end














return;