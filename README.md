# SpikeTools for Takahashi Lab
@ Graduate School of Brain Science, Doshisha University

***Tools for analyzing neural activity recorded using our custom-made software.***
***The software was developed using intan-tech libraries for LabVIEW programming language.***

## MATLAB path setup

Most legacy top-level MATLAB functions have been moved into `legacy/` to keep
the repository root clean. Before using the toolbox, run this once per MATLAB
session from the repository root:

```matlab
setup_spiketools_path
```

This adds the repository and all subfolders to the MATLAB path, so older calls
such as `preProcessIlvrc`, `batchSpikesIntan`, and `extractPos` continue to work.

## Repository layout

Top-level folders are organized by analysis domain or acquisition system:

| Folder | Contents |
|---|---|
| `analysis/` | Shared analysis helpers, including LFP and statistics utilities. |
| `configFiles/` | Configuration files and electrode/tetrode layout files. |
| `event/` | Event extraction utilities for noisy event channels and RHX data. |
| `extractSpk/` | Spike extraction MEX files and related source code. |
| `fish/` | Fish behavior and magnetic-field experiment analysis pipeline. |
| `HeadDirection/` | Head-direction tuning, polar map, occupancy, and circular statistics code. |
| `LFP/` | LFP filtering, wavelet, phase, and raster-related functions. |
| `mouse/` | Mouse-specific preprocessing, RF/event extraction, and LFP analysis code. |
| `neuropixels/` | Neuropixels probe bank configuration and IMRO helpers. |
| `PD/` | Photodiode/visual-stimulus response analysis routines. |
| `Place/` | Place-cell, place-map, occupancy, lap, and spatial information analysis. |
| `procedure/` | Procedure-level scripts for specific animal or experiment workflows. |
| `SGLX/` | SpikeGLX-related spike matrix and conversion helpers. |
| `spikegadgets/` | SpikeGadgets/Trodes export notes, notebooks, and Phy conversion helpers. |
| `spikesorting/` | KiloSort, MountainSort, cluster quality, waveform, and sorting conversion tools. |
| `Time/` | Time-cell and sequence-map analysis functions. |
| `turtle/` | Sea turtle batch scripts and result summary scripts. |

`legacy/` contains older MATLAB entry points that used to live in the repository
root. They are still usable after running `setup_spiketools_path`.

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
