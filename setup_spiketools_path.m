function setup_spiketools_path()
%SETUP_SPIKETOOLS_PATH Add this repository and all subfolders to MATLAB path.
repoRoot = fileparts(mfilename('fullpath'));
paths = regexp(genpath(repoRoot), pathsep, 'split');
paths = paths(~cellfun(@isempty, paths));
keep = cellfun(@(p) isempty(strfind(p, [filesep '.git'])), paths);
addpath(strjoin(paths(keep), pathsep));
fprintf('SpikeTools path added: %s\n', repoRoot);
end
