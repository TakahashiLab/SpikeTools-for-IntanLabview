function [cOut, kkOut] = coinSpk(kkOut, msT, fstart, Traj, varargin)
    p = inputParser;
    p.addParamValue('timeWindow', 1, @isnumeric);
    p.addParamValue('samplingRate', 31.25, @isnumeric);
    p.addParamValue('target', -1, @isnumeric);

    p.parse(varargin{:});

    % 同時スパイクを定義する時間ウィンドウ（ミリ秒）
    timeWindow = p.Results.timeWindow; % ±5ms
    samplingRate = p.Results.samplingRate * 10 ^ 3; % 31.25kHz
    kHz = p.Results.samplingRate;
    timeWindowSamples = timeWindow * 1e-3 * samplingRate; % サンプル単位での時間ウィンドウ
    targetNeuron = p.Results.target; %target neuron #

    if size(Traj, 1) ~= size(msT, 1)
        msFPS = median(diff(msT));
        seq = 1:length(Traj);
        seq = seq - 1;
        msT = seq * msFPS + msT(1);
        msT = msT';
    end

    lastKK = cell(1, 3);
    lastKK{1, 3} = ceil((msT' - fstart) .* kHz);
    lastKK{1, 1} = zeros(4, length(lastKK{1, 3}));
    kkOut = [kkOut; lastKK];

    step = size(kkOut{1, 1}, 2) / length(kkOut{1, 3});

    % 各ニューロンのペアに対して処理
    numNeurons = size(kkOut, 1);
    cOut = cell(1, 5); % 同時スパイクを格納するセル配列
    oOut = cOut;
    counter = 1;
    ts = []; %spikes of targetNeuron fired together with others

    % 全てのニューロンのスパイクタイムスタンプを集約
    allSpikes = cell(numNeurons, 1);

    for n = 1:numNeurons
        allSpikes{n} = kkOut{n, 3};
    end

    % 二つのニューロン間の全てのスパイクタイムスタンプの組み合わせを計算
    for i = 1:numNeurons

        for j = i + 1:numNeurons
            spikes1 = allSpikes{i};
            spikes2 = allSpikes{j};

            % 各スパイクペア間の時間差を計算
            timeDiffs = abs(bsxfun(@minus, spikes1, spikes2'));

            % 同時スパイクを識別
            [row, col] = find(timeDiffs <= timeWindowSamples);

            for k = 1:numel(row)
                s1 = col(k);
                s2 = row(k);
                % 同時スパイクをcOutに格納

                cOut{counter, 1} = [cOut{counter, 1} getSpikesD(kkOut{i, 1}, step, s1)]; % ニューロンiの波形
                cOut{counter, 2} = [cOut{counter, 2} getSpikesD(kkOut{j, 1}, step, s2)]; % ニューロンjの波形
                cOut{counter, 3} = [cOut{counter, 3} spikes1(s1)]; % スパイクタイミング
                cOut{counter, 4} = [cOut{counter, 4} spikes2(s2)];
                cOut{counter, 5} = [i j]; % ニューロンのペア

                if targetNeuron

                    if any([i j] == targetNeuron)

                        if i == targetNeuron
                            ts = [ts s1];
                        elseif j == targetNeuron
                            ts = [ts s2];
                        end

                    end

                end

            end

            if ~isempty(row)
                counter = counter + 1;
                cOut = [cOut; oOut];
            end

        end

    end

    if isempty(cOut{counter, 1})
        cOut(counter, :) = [];
    end

    if targetNeuron

        if ~isempty(ts)
            singleSpikes = setdiff(1:length(kkOut{targetNeuron, 3}), ts);
            kkOut{targetNeuron, 3} = kkOut{targetNeuron, 3}(singleSpikes);
            kkOut{targetNeuron, 1} = getSpikesD(kkOut{targetNeuron, 1}, step, singleSpikes);
        end

    end
