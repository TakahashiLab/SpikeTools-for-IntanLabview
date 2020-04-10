# SpikeTools

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
  
# How to use

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
