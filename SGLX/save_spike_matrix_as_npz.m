function save_spike_matrix_as_npz(Out,filename_prefix)
    % データ準備
   

    % スパイクタイミング（Outの3列目にスパイクタイミングが格納されている）
    spike_times_cell = Out(:, 3);  % セル配列としてスパイクタイミングを抽出

    % パラメータ設定
    bin_width = 6000;  % 200msを秒単位に変換, 30kHz
    start_time = 0;    % 計測開始時間
    end_time = 0;      % 終了時間の初期化

    % ニューロン数を取得
    num_neurons = length(spike_times_cell);

    % 各ニューロンのスパイクタイミングから最大時間を取得
    for i = 1:num_neurons
        if ~isempty(spike_times_cell{i})  % スパイクタイミングが存在する場合
            end_time = max([end_time, max(spike_times_cell{i})]);  % 最大時間を更新
        end
    end
    
    % 時間ビンの数を計算
    num_bins = ceil((end_time - start_time) / bin_width);

    % スパースマトリックスを初期化
    spike_matrix_sparse = sparse(num_neurons, num_bins);

    % 各ニューロンごとにスパイクタイミングを処理
    for i = 1:num_neurons
        if ~isempty(spike_times_cell{i})  % スパイクタイミングが存在する場合
            spike_times = double(spike_times_cell{i});  % スパイクタイミングを数値に変換
            bin_indices = ceil(spike_times / bin_width);  % スパイクタイミングをビンに変換
            spike_matrix_sparse(i, bin_indices) = spike_matrix_sparse(i, bin_indices) + 1;  % スパイクマトリックスに記録
        end
    end

    % スパースマトリックスをフルマトリックスに変換
    spike_matrix_full = full(spike_matrix_sparse);
    


    % ファイル名を作成
    output_filename = sprintf('%s.npz', filename_prefix);

    % .npzファイルとして保存 (Pythonとの連携が必要)
      py.numpy.savez(output_filename, spks=spike_matrix_full);

    disp(['Spike matrix saved as ', output_filename]);
end
