% on the target neuropixels data diretory
% [Out,msT,st,Pos,Sts,vt]=readSpks('.','animal','turtle');
% save behav.mat

function results=batchSeaturtleGPT(Out, msT, Pos, Sts,varargin)
p = inputParser;
p.addParamValue('animal', 'rat', @ischar);
p.addParamValue('shufflen', 100, @isnumeric);
p.addParamValue('shuffletype', 'shift10m', @ischar);
p.addParamValue('folded', 4, @isnumeric);
p.parse(varargin{:});
animal = p.Results.animal;
shuffleN = p.Results.shufflen;
shuffleType = p.Results.shuffletype;
folded = p.Results.folded;%for seaturtle, 4


bothside=0;
% 条件設定
condN=size(Sts,2);
if condN==4
    conditions = {Sts{1}, Sts{2}, Sts{3}, Sts{4}};
elseif condN==2
    conditions = {Sts{1}, Sts{2}};
end
num_neurons = size(Out, 1);

% 結果を格納する構造体
results = struct();
results.individual = zeros(num_neurons, condN*3); % 各条件での有意性
results.normal_field = zeros(num_neurons, 4); % 通常磁場での有意性
results.zero_field = zeros(num_neurons, 4); % ゼロ磁場での有意性
results.peakfr = zeros(num_neurons,condN);
retuls.ratemap = zeros(num_neurons,condN,60);

% 各ニューロンについて処理
for neuron_idx = 1:num_neurons
    fprintf('%d/%d\n',neuron_idx,num_neurons);
    spikes = Out{neuron_idx, 3}; % 例としてスパイクタイミングが3列目に格納されていると仮定
    for cond_idx = 1:condN
        condition = conditions{cond_idx};
        pos = Pos(condition, :);
        time = msT(condition);

        % makeSigShuffleを呼び出して計算
        [p0, ps, theta] = makeSigShuffle(spikes, pos, time, 'animal', 'seaturtle', 'shufflen', 100,'shuffleType',shuffleType,'folded',folded);
        mr0=p0{2};
        mvl0=p0{3};
        ratemap=p0{1};
        %bidirectional
        mr0b=p0{6};
        mvl0b=p0{7};

        mrs=ps{2};
        mvls=ps{3};

        %bidirectional
        mrsb=ps{5};
        mvlsb=ps{6};

        ang_hd0=p0{4};
        peakfr0=p0{5};

        % 有意性の検定
        mrs=sort(mrs);
        mvls=sort(mvls);
        mr0_significance = mr0 > mrs(95);
        mvl0_significance = mvl0 > mvls(95);
        %bidirectional
        mrsb=sort(mrsb);
        mvlsb=sort(mvlsb);
        mr0b_significance = mr0b > mrsb(95);
        mvl0b_significance = mvl0b > mvlsb(95);

        results.individual(neuron_idx, cond_idx) = mr0_significance;
        results.individual(neuron_idx, cond_idx+condN) = mvl0_significance;
        results.individual(neuron_idx, cond_idx+condN*2) = ang_hd0;
        results.individual(neuron_idx, cond_idx+condN*3) = mr0b_significance;
        results.individual(neuron_idx, cond_idx+condN*4) = mvl0b_significance;
        results.peakfr(neuron_idx,cond_idx) = peakfr0;
        results.ratemap(neuron_idx,cond_idx,:)=ratemap;

        if condN==4
            % 通常磁場の有意性
            normal_conditions = [Sts{1}; Sts{3}];
            pos_normal = Pos(normal_conditions, :);
            time_normal = msT(normal_conditions);

            [p0, ps, theta] = makeSigShuffle(spikes, pos_normal, time_normal, 'animal', 'seaturtle', 'shufflen', 100,'shuffleType',shuffleType,'folded',folded);
            mr0=p0{2};
            mrs=ps{2};
            ang_hd0=p0{4};
            peakfr0=p0{5};
            ratemap=p0{1};

            mrs=sort(mrs);

            %bidirectional
            mr0b=p0{6};
            mvl0b=p0{7};

            mvl0=p0{3};
            mvls=ps{3};
            mvls=sort(mvls);
            %bidirectional
            mrsb=ps{5};
            mvlsb=ps{6};

            mrsb=sort(mrsb);
            mvlsb=sort(mvlsb);

            mr0_significance_normal = mr0 > mrs(95);
            mvl0_significance_normal = mvl0 > mvls(95);
            mr0b_significance_normal = mr0b > mrsb(95);
            mvl0b_significance_normal = mvl0b > mvlsb(95);

            results.normal_field(neuron_idx, 1) = mr0_significance_normal;
            results.normal_field(neuron_idx, 2) = mvl0_significance_normal;
            results.normal_field(neuron_idx, 3) = ang_hd0;
            results.normal_field(neuron_idx, 4) = peakfr0;
            results.normal_field(neuron_idx, 5) = mr0b_significance_normal;
            results.normal_field(neuron_idx, 6) = mvl0b_significance_normal;

            results.ratemap(neuron_idx,4+1,:)=ratemap;

            % ゼロ磁場の有意性
            zero_conditions = [Sts{2}; Sts{4}];
            pos_zero = Pos(zero_conditions, :);
            time_zero = msT(zero_conditions);

            [p0, ps, theta] = makeSigShuffle(spikes, pos_zero, time_zero, 'animal', 'seaturtle', 'shufflen', 100,'shuffleType',shuffleType,'folded',folded);
            mr0=p0{2};
            mrs=ps{2};
            ang_hd0=p0{4};
            peakfr0=p0{5};
            ratemap=p0{1};

            mrs=sort(mrs);
            mvl0=p0{3};
            mvls=ps{3};
            mvls=sort(mvls);
            %bidirectional
            mr0b=p0{6};
            mvl0b=p0{7};
            %bidirectional
            mrsb=ps{5};
            mvlsb=ps{6};

            mrsb=sort(mrsb);
            mvlsb=sort(mvlsb);

            mr0_significance_zero = mr0 > mrs(95);
            mvl0_significance_zero = mvl0 > mvls(95);
            mr0b_significance_zero = mr0b > mrsb(95);
            mvl0b_significance_zero = mvl0b > mvlsb(95);
            results.zero_field(neuron_idx,1) = mr0_significance_zero;
            results.zero_field(neuron_idx,2) = mvl0_significance_zero;
            results.zero_field(neuron_idx,3) = ang_hd0;
            results.zero_field(neuron_idx,4) = peakfr0;
            results.zero_field(neuron_idx,5) = mr0b_significance_zero;
            results.zero_field(neuron_idx,6) = mvl0b_significance_zero;

            results.ratemap(neuron_idx,4+2,:)=ratemap;
        end
    end


end
