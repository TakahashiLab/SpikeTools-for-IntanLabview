function [rensemble] = rawCellClass(rensemble, varargin)
    p = inputParser;
    p.addParamValue('raw', [], @ismatrix);
    p.addParamValue('step', 32, @isnumeric);
    p.addParamValue('samplingrate', 25000, @isnumeric);

    p.parse(varargin{:});
    x = p.Results.raw;
    step = p.Results.step;
    SamplingRate = p.Results.samplingrate;
    OneMs = round(SamplingRate / 1000);

    for i = 1:size(rensemble, 1)
        fprintf('Unit %d\n', i);
        bbwf = rensemble{i, 1};
       
        mbbwf = reshape(bbwf, size(bbwf, 1), step, size(bbwf, 2) / step);
        mbbwf = permute(mbbwf, [3 1 2]);

        mbbwf = squeeze(mean(mbbwf)) - mean(mean(mean(mbbwf)));

        for t = 1:size(mbbwf, 1)
            [a(t), b(t)] = max(abs(mbbwf(t, :)));
        end

        [aa, bb] = max(a, [], 2);
        MaxWaves = mbbwf(bb, :);
        %% get trough-peak delay times
        %AllWaves(:,:,1) = [];

        thiswave = MaxWaves';
        [minval, minpos] = min(thiswave);
        minpos = minpos(1);

        [maxval, maxpos] = max(thiswave);
        [dummy, maxpos] = max(thiswave(minpos + 1:end));

        if isempty(maxpos)
            warning('Your Waveform may be erroneous')
            maxpos = 1;

        end

        maxpos = maxpos(1);
        maxpos = maxpos + minpos;
        tp = maxpos - minpos; %In number of samples
        tp = tp / OneMs;
        rensemble{i, 4} = tp;

        clear aa bb a b;

        %% get spike width by taking inverse of max frequency in spectrum (based on Adrien's use of Eran's getWavelet)

        w = thiswave;
        w = [w(1) * ones(1000, 1); w; w(end) * ones(1000, 1)];
        [wave f t] = getWavelet(w, SamplingRate, 500, 3000, 128);
        %We consider only the central portion of the wavelet because we
        %haven't filtered it before hand (e.g. with a Hanning window)
        wave = wave(:, int16(length(t) / 4):3 * int16(length(t) / 4));
        %Where is the max frequency?
        [maxPow ix] = max(wave);
        [dumy mix] = max(maxPow);
        ix = ix(mix);
        spkW = 1000 / f(ix);

        
        rensemble{i, 5} = spkW;
    end

end
