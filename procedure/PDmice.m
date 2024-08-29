%set working directory
fn=pwd;
%
%set reference
%split hemisphere
forceRef{1,1}=1;forceRef{1,2}=1:4;forceRef{2,2}=[5:8];
%%
%load images....mat
Traj=FixPos(Pos);
save DLC.mat Traj;
%batch process for constructing neural and event signals
%1fn:working directory
%2: process option, 0:all
%3: reference assignment: forceRef
%4: verion: 2018 or new
%5: LED position: tetrode number
%6: realtime feedback?: 
%7  gpu to use?: 
%example for continuous stimulation
preProcessIlvrc4Mouse(fn,0,forceRef,0,3,0,0);
%example for realtime feedback
preProcessIlvrc4Mouse(fn,0,forceRef,0,3,1,0);

%%%%%%%%%%%%%%%
%%make ensemble file
%load results from klustakwik
fn=pwd;
kkOuts=batchLoadKK(fn);
%make ensemble file
%for example, choose clusters #1-4 and the merge #1,2 and #3,4
an{8}=[];
en{8}=[];
an{1}=[1 2 3 4];
en{1}={[1 2];[3 4]};
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


%cell identification
load tmp.mat;
[rkkOuts,rensemble]=makeBroadBand(fn,'raw',x);
rensemble=rawCellClass(rensemble);
save rensemble.mat rensemble rkkOuts;
[phaseHistPyr,phaseHistInt,PyrIntList,PyrIntListStim,fr,tp,sw,pi,phaseHistPyrCtrl, phaseHistIntCtrl,cq]=contPD('2018','TK18052901','method','cellclassify');


%cell classification
%
%This program produces cellclass.mat with cell classification data.
batchCondPD('deepMachine','cellclassify');
[TPs,SWs,PIs]=Cellall(fn);

%move cellclass.mat to mouse/Gain
%then cd mouse/Gain
%
batchCondPD('deepMachine','gainmap');% deepMachine=ubuntu, windows=windows server

%Gainmap
[phaseHistPyr,phaseHistInt,PyrIntList,PyrIntListStim,fr,tp,sw,pi,phaseHistPyrCtrl, phaseHistIntCtrl]=contPD('2019','TK19050801','cellClass',0,'localcell',1);

%led or correlation tagging, ascending wave, gamma-band, neuron by phase plot, histogram, top frequency 120Hz, cell type cc(cross-correlation / opto-tagging) 
Gainall(fn,'stimulation','ledtagstim','display','ascend','plot','all','band','gamma','xyaxis','neuronphase','hist','on','topFreq',120,'celltype','cc');

%Coherence
Cohereall(fn,'stimulation','ledtagstim','display','both','plot','all','band','all','xyaxis','freqphase','hist','on','topFreq',40,'celltype','all');