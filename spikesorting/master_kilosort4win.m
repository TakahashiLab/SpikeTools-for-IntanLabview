function master_kilosort4win(filename)

if nargin==0
    filename='test2.bin';
end

% default options are in parenthesis after the comment
homefn='c:/users/stakahas/Documents/Matlab/';% change this directory first
KSfn=fullfile(homefn,'KiloSort');
addpath(genpath(KSfn)); % path to kilosort folder
NPYfn=fullfile(homefn,'npy-matlab');
addpath(genpath(NPYfn)); % path to npy-matlab scripts

%pathToYourConfigFile = '/home/tsusumu/KiloSort/configFiles'; % take from Github folder and put it somewhere else (together with the master_file)
pathToYourConfigFile = fullfile(homefn,'SpikeTools-for-IntanLabview/configFiles');

%from StandardConfig_MOVEME.m
cwfn=pwd;
ops.fbinary             = fullfile(cwfn,filename); % will be created for 'openEphys'		
ops.fproc               = fullfile(cwfn,'temp_wh.dat'); % residual from RAM of preprocessed data		
ops.root                = cwfn;
ops.chanMap             = fullfile(pathToYourConfigFile,'tetrode64.mat'); % make this file using createChannelMapFile.m	

run(fullfile(pathToYourConfigFile, 'StandardConfig_MOVEME.m'))

tic; % start timer
%
if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

if strcmp(ops.datatype , 'openEphys')
   ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
end
%
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and
                                          % extract spikes for
                                          % initialization

rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% AutoMerge. rez2Phy will use for clusters the new 5th column of st3 if you run this)
%     rez = merge_posthoc2(rez);

% save matlab results file
save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');

% save python results file for Phy
rezToPhy(rez, ops.root);

% remove temporary file
delete(ops.fproc);
%% 
%