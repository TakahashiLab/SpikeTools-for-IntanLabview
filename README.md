# SpikeTools for Takahashi Lab
@ Graduate School of Brain Science, Doshisha University

***Tools for analyzing neural activity recorded using our custom-made software.***
***The software was developed using intan-tech libraries for LabVIEW programming language.***

- Please clone the following tools from cortexlab at UCL to your home directory (c:/users/[your user account here]/Documents/Matlab) using git. 
  * [KiloSort](https://github.com/cortex-lab/KiloSort)
  * [spikes](https://github.com/cortex-lab/spikes)
  * [npy-matlab](https://github.com/kwikteam/npy-matlab)
  * [chronux](https://github.com/jsiegle/chronux)

- Copy compiled files in the KiloSortCUDAfolder at TakahashiLabB HDD into the respective directory of your cloned KiloSort directory.  

- Copy startup.m file in the MATLABpathfolder at TakahashiLabB HDD into your MATLAB home directory.

- Curation 
  * [phy](https://github.com/cortex-lab/phy)
---  

# How to use with KlustaKwik or ICSort and tTools  

1. Convert your ilvrc file into MAT files.
```matlab
fn=pwd;
preProcessIlvrc(fn,0);
```

1'. PD mouse
```matlab
fn=pwd;
load DLCAnalyzedTrajectory.mat
Traj=Pos;
save DLC.mat;
forceRef{1,1}=1;forceRef{1,2}=1:4;forceRef{2,2}=[5:8];%split median referencing
preProcessIlvrc4Mouse(fn,0,forceRef,0,LEDp,0);% 0:entire process, 0:old version/ 1:new version, LEDp:LED position (tetrode #), realtime feedback? No/Yes (0/1) 
```

2. Spike sorting using KlustaKwik.
```matlab
fn=pwd;
batchSpikesIntan(fn);
```

3. Extract the timing of video frames from event  
```matlab
fn=pwd;
PosT = extractPos(fn,Traj);
```

4. Construct place map of unit #3
```matlab
[rate_map, ~, ~, oc_map] = pmap(kkOut{3,3},Traj,PosT,0,'animal','rat');
imagePmap(rate_map,oc_map);
```

5. Construct head direction polar map of unit #3
```matlab
HeadDirectionPlot = plot_polar_rate_map(kkOut{3,3},Traj,PosT,0,'rat');
```

---
# How to use with KiloSort and Phy

1. Convert your ilvrc file into .bin file.
```matlab
fn=pwd;
preProcessIlvrc(fn,4);
```

2. Spike sorting using KiloSort.
```matlab
master_kilosort4win('hogehoge.bin');
```

3. Curation with phy on Anaconda-prompt
```matlab
phy template-gui ~/hogehoge/hoge/params.py
```

4. Convert KiloSort timestamps into sampling counts
```matlab
fn=pwd;
ks=convertKS(fn);
```

5. Extract the timing of video frames from event  
```matlab
fn=pwd;
PosT = extractPos(fn,Traj);
```

6. Construct place map of unit #3
```matlab
[rate_map, ~, ~, oc_map] = pmap(ks.st(ks.clu==3),Traj,PosT,0,'animal','rat');
imagePmap(rate_map,oc_map);
```

7. Construct head direction polar map of unit #3
```matlab
HeadDirectionPlot = plot_polar_rate_map(ks.st(ks.clu==3),Traj,PosT,0,'rat');
```

---
# How to use with mountainlab (mountainsort and mountainview) 

1. Convert your e.mat file into .mda file.
```matlab
fn=pwd;
mountainWrite(fn,16);%for 16 tetrodes
```

2. Spike sorting using mountainsort.
```python
conda activate mountain
cd 'data diretory'
python /nfs/home/tsusumu/SpikeTools-for-IntanLabview/spikesorting/master_ms4.py
```

3. Curation using qt-mountainview
```python
python /nfs/home/tsusumu/SpikeTools-for-IntanLabview/spikesorting/master_mv.py
```

4. Convert mountainsort timestamps into sampling counts
```python
fn=pwd;
[ensemble,an,en]=convertMS(fn,16);%for 16 tetrodes
```
