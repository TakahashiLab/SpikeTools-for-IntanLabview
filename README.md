# SpikeTools
Tools for analyzing neural activity recorded using our custom-made software.
The software was developed using intan-tech libraries for LabVIEW programming language.

1. Please clone the following tools from cortexlab at UCL to your home directory (c:/users/[your user account here]/Documents/Matlab) using git. 
* KiloSort (https://github.com/cortex-lab/KiloSort)
* phy (https://github.com/cortex-lab/phy)
* spikes (https://github.com/cortex-lab/spikes)
* npy-matlab (https://github.com/kwikteam/npy-matlab)
* chronux (https://github.com/jsiegle/chronux)

2. Copy compiled files in the KiloSortCUDAfolder at TakahashiLabB HDD into the respective directory of your cloned KiloSort directory.  

3. Copy startup.m file in the MATLABpathfolder at TakahashiLabB HDD into your MATLAB home directory.

[How to use]

Convert your ilvrc file into .bin file.
> fn=pwd;preProcessIlvrc(fn,4);

Spike sorting using KiloSort.
> master_kilosort4win('hogehoge.bin');

Curation with phy
on Anaconda-prompt
> phy template-gui ~/hogehoge/hoge/params.py

Convert KiloSort timestamps into sampling counts
> fn=pwd;ks=convertKS(fn);

Extract the timing of video frames from event  
> fn=pwd; PosT = extractPos(fn,Traj);

Construct place map of unit #3
> [rate_map, ~, ~, oc_map] = pmap(ks.st(ks.clu==3),Traj,PosT,0,'animal','rat');
> imagePmap(rate_map,oc_map);

Construct head direction polar map of unit #3
> HeadDirectionPlot = plot_polar_rate_map(ks.st(ks.clu==3),Traj,PosT,0,'rat');

