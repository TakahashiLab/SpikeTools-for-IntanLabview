function [rkkOuts,rensemble,IS,LR] = makeBroadBand(basename, varargin)
    p = inputParser;
    p.addParamValue('raw', [], @ismatrix);
    p.addParamValue('step', 32, @isnumeric);

    p.parse(varargin{:});
    x = p.Results.raw;
    step = p.Results.step;

    %load all data
    [path, name, ext] = fileparts(basename);

    dataFolder = fullfile(path, name);

    if isempty(x)
        fn = fullfile(dataFolder, ['tmp.mat']);
        load(fn, 'x');
    end

    fn = fullfile(dataFolder, ['ensemble.mat']);
    load(fn, 'an', 'en', 'kkOuts');
    numOfElectrodes = size(x, 1) / 4;

    for i = 1:numOfElectrodes
        fprintf('Processing tetrode number %d\n', i);
        bbwf = double(x(1 + 4 * (i - 1):4 * i, :)); %broad band tetrode waveform

        for j = 1:length(an{i})
            cn = an{i}(j); %cluster number
            rkkOuts{i}{cn, 3} = kkOuts{i}{cn, 3};
            %added modified for first spikes
            rkkOuts{i}{cn, 3}(kkOuts{i}{cn, 3} < step) = [];
            
            rkkOuts{i}{cn, 1} = extractMS(bbwf, step, double(rkkOuts{i}{cn, 3}));

        end

    end
    [rensemble,IS,LR]=makeEnsemble64(rkkOuts,an,en);
    rensemble=rawCellClass(rensemble);
    
end
