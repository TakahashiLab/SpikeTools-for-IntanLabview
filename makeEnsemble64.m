function [ensemble, IS, LR] = makeEnsemble64(kkOuts, an, en, varargin)
    p = inputParser;
    p.addParamValue('step', 32, @isnumeric);
    p.addParamValue('khz', 25, @isnumeric);

    p.parse(varargin{:});
    step = p.Results.step;
    kHz = p.Results.khz;

    ensemble = [];
    IS = [];
    LR = [];
    tetNum = size(kkOuts, 2);

    for i = 1:tetNum
        tmp = kkOuts{i}(an{i}, :);

        if ~isempty(tmp)
            step = size(tmp{1, 1}, 2) / size(tmp{1, 3}, 2);
        end

        if ~isempty(en{i})

            rInd = [];

            for j = 1:size(en{i}, 1)
                sameCell = sort(en{i}{j});

                oind = find(an{i} == sameCell(1));
                n = tmp(oind, :);

                for k = 2:length(sameCell)
                    ind = find(an{i} == sameCell(k));
                    rInd = [rInd ind];
                    n{1, 1} = [n{1, 1} tmp{ind, 1}];
                    n{1, 3} = [n{1, 3} tmp{ind, 3}];

                    %reorder spike timings and waveform
                    [n{1, 3}, ind] = sort(n{1, 3});
                    n{1, 1} = getSpikesD(n{1, 1}, step, ind);
                end

                %            tmp=[tmp;n];
                tmp(oind, :) = n;
            end

            tmp(rInd, :) = [];
        end

        ensemble = [ensemble; tmp];
        fprintf('tetNum=%d\n', i);

        if size(tmp, 1) == 1
            is = NaN;
            lr = NaN;
        else
            [is, lr] = clusterQualitys(tmp, step, kHz);
        end

        IS = [IS; is(:, 1)];
        LR = [LR; lr(:, 1)];
    end

    return;
