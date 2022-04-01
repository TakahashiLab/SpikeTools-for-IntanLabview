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
%tetrode #1
%
j=1;
for i=1:size(kkOuts{j},1)
    figure;
    dispDodeca(kkOuts{j},i,'limitlen',100);
end

ensemble=makeEnsemble64(kkOuts,an,en);
%Please save ensemble information for later analyses
save ensemble.mat ensemble an en kkOuts;
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
[NosePork,Treadmill,Pos,PosT]=extractRHX(event,t,Pos);   

save event.mat NosePork Treadmill PosT;
load DLC.mat
Pos=Traj;
save positions.mat Pos PosT;
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
load ensemble.mat ensemble
for i=1:size(ensemble,1)
    figure;
    % reference point: nose-porking #2 
    % pre- and post duration: 8 sec
    % smoothing on
    plotTime(ensemble{i,3},beNP{2},'verbose',1,'pre',3,'post',7,'smooth', ...
             1);
end
%%%%%%%%%%%%%%
[beNP,tmr,tmSpeed]=sawataniR1(NosePork,Treadmill);
%%%%%%%%
[seq,order]=SequenceTmap4S(ensemble,beNP{3},'jitterpre',10, ...
                           'jitterpost',7);
imagesc(seq(order,:));