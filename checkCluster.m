function checkCluster(ensemble,an,en)
tm=tetrodeMap(an,en);
demoFeatureCore(ensemble,tm,tetNum);
return;
%%%%%%%%%%%%%%%%
function demoFeatureCore(ensemble,tm,tetNum)

kHz=30;
step=size(ensemble{1,1},2)/size(ensemble{1,3},2);
c={'r.','g.','b.','k.','y.','m.','r*','g*','k*','y*'};

FetName={'PC1','PC2','PC3','Energy'};
tetN=size(ensemble{1,1},1);
FetN=tetN*(3+1);
N=4;

    tInd=find(tm==tetNum);
    
    
    if length(tInd)>1
        units=ensemble(tInd,:);
        %[is,lr]=clusterQualitys(units,step,kHz);
        %is=is(:,1);
    elseif length(tInd)==1
        units=ensemble(tInd,:);
    end
    
        loop=size(units,1);
        for i=1:loop
            wf=units{i,1};
            Fets{i}=getFets(wf,step);
        end    

        fprintf('tetrode#%d\n',tetNum);
        for i=1:(FetN-1)
            figure;
            cn=1;
            for j=(i+1):FetN
                subplot(4,4,cn)
                cn=cn+1;
                tn1=ceil(i/4);
                fn1=mod(i,4);
                fn1(fn1==0)=4;
                tn2=ceil(j/4);
                fn2=mod(j,4);
                fn2(fn2==0)=4;
                tn=sprintf('t%d-%s vs t%d-%s',tn1,FetName{fn1},tn2,FetName{fn2});
                hold on;
                title(tn);
                for k=1:loop
                    plot(Fets{k}(:,i),Fets{k}(:,j),c{k});
                end
            end
        end

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fets=getFets(wf,step)
ewf=sqrt(wf.^2);
Energy=[];
lenTs=size(wf,2)/step;


for i=1:lenTs
  tmpEnergy=sum(ewf(:,1+(i-1)*step:i*step),2)/step;
  Energy=[Energy tmpEnergy];
  tmpEnergy(find(tmpEnergy==0))=1;
end



PCA=makePCA(wf,size(wf,1),step);

Fets=[PCA Energy'];

return;