%save  Deeplabcut-processedHdf5.mat as DLC.mat
%%%%
fn=pwd;
forceRef{1,1}=0;
preProcessRHX(fn,0,forceRef,1);
batchSpikesIntan(fn,5,1:16,50,'rhd');
kkOuts=batchLoadKK(fn);
%make ensemble file
%for example, choose clusters #1-4 and the merge #1,2 and #3,4
an{16}=[];
en{16}=[];
an{1}=[1 2 3 4];
en{1}={[1 2];[3 4]};
%%%%
%%%%% check spike waveform
%%%
fn=pwd;
kkOuts=batchLoadKK(fn);
%
%check tetrode #16 with respect to spike waveform.
%
wfCheck(kkOuts,16);

ensemble=makeEnsemble64(kkOuts,an,en);
%Please save ensemble information for later analyses
save ensemble.mat ensemble an en kkOuts -v7.3;

%
%find neurons from ensemble
%for example,you can find neurons in the tetrode #1 
tet1=ensemble(find(tetrodeMap(an,en)==1),:);
%
%%%%%%%%%%%%%%%%%%%%%%%%
%%check cluster in the feature space
%ex. check cluster of tetrode #1 in the feature space 
checkClusterM(ensemble,an,en,1)

load ???Event.mat

    load DLC.mat Pos;
    Pos=Pos';
    [NosePork,Treadmill,Pos,PosT]=extractRHX(event,t,Pos);   
    PosT=PosT';

save positions.mat Pos PosT;
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% plot Time field
load ensemble.mat ensemble
for i=1:size(ensemble,1)
    figure;
    % reference point: nose-porking #2 
    % pre- and post duration: 8 sec
    % smoothing on
    plotTime(ensemble{i,3},beNP{2},'verbose',1,'jitterpre',3,'jitterpost',7,'smooth', ...
             1);
end
%%%%%%%%%%%%%% plot theta phase precession
tetrodePoint=tetrodeMap(an,en);
samplrate=1000;
for i=1:size(ensemble,1)
    figure;
    ts=thetaPhase(dlfp(1+(tetrodePoint(i)-1)*4,:),samplrate);    
    plotTimePhase(ensemble{i,3},beNP{2},ts,'verbose',1','jitterpre',2,'jitterpost',10,'smooth',0);
end

%%%%%%%%%%%%%%
[beNP,tmr,tmSpeed]=sawataniR1(NosePork,Treadmill);
%%%%%%%%sequential activity pattern on time dimension
[seq,order]=SequenceTmap4S(ensemble,beNP{3},'jitterpre',10, ...
                           'jitterpost',7);
imagesc(seq(order,:));
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
%%%%%% linearlized route
[LTraj,TMpoints,GI]=fitLinearModel(Pos,'posture',[5 6],'mazetype', ...
                                   'linear');
%%%%%%%%%
%%%%%% exclude treadmill (delay) period
[ensemble2,LTraj2,PosT2]=excludeDelay(ensemble,beNP{2},beNP{3}, ...
                                      LTraj,PosT);
%%%%%%%%%%
%%%%%%% sequential activity pattern on space dimension
[seq,order,oc_map]=SequencePmap4S(ensemble2,LTraj2,PosT2);
imagesc(seq(order,:));

%%%%%%%%%%%%%%%%%%%%%timepoints of forward/backward laps
[fTP,bTP]=extractLap(PosT,beNP{1},beNP{3});
%plot(LTraj);
%hold on;
%plot(fTP,LTraj(fTP),'r.');
%plot(bTP,LTraj(bTP),'g.');