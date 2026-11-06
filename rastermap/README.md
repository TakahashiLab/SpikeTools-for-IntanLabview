# Rastermap for curated good clusters

Phyで`good`に分類したKilosortクラスタだけを集め、Rastermap用の
`neurons × time bins`行列を作るスクリプトです。1回の録画に含まれる
複数プローブ／シャンクの`imec*_ks4`をまとめて処理できます。

## 入力

指定したルート以下から、次を含む`imec*_ks4`フォルダを再帰検索します。

- `spike_times.npy`
- `spike_clusters.npy`
- `cluster_group.tsv`
- `params.py`（なくても既定の30 kHzで処理）

`cluster_group.tsv`の`group`、`KSLabel`または`label`列が`good`の
クラスタだけを使用します。バックアップ用の`*_orig`フォルダは除外します。

## インストール

PowerShellまたはAnaconda Promptで実行します。

```powershell
cd "J:\OneDrive - 同志社大学\Documents\MATLAB\SpikeTools-for-IntanLabview"
python -m pip install -r rastermap\requirements.txt
```

Rastermapを実行せず行列作成だけ行う場合は、`numpy`と`matplotlib`だけでも
動作します。

## 標準的な使い方

各録画の`*_imec0`フォルダへ移動して、ルートに`.`を指定します。
以下は既存データと同じ100 ms bin、平方根変換後z-scoreの設定です。

```powershell
cd "P:\seabird\2025\original\subjectW\20250918\catgt\catgt_2025091801J_g0\2025091801J_g0_imec0"

python "J:\OneDrive - 同志社大学\Documents\MATLAB\SpikeTools-for-IntanLabview\rastermap\run_good_clusters_rastermap.py" `
  . `
  --out .\rastermap `
  --bin-ms 100 `
  --run-rastermap
```

この例では、録画フォルダ内の4つの`imec*_ks4`をまとめて処理します。
Rastermapをまだインストールしていない場合、まず`--run-rastermap`を
外せば入力変換と未整列プレビューだけを作成できます。

複数のルートを明示することもできます。

```powershell
python .\rastermap_good_clusters\run_good_clusters_rastermap.py `
  "P:\path\to\recording_A_imec0" `
  "P:\path\to\recording_B_imec0" `
  --out "P:\seabird\rastermap_combined" `
  --bin-ms 100 `
  --min-spikes 100 `
  --run-rastermap
```

注意：異なる録画を同時指定した場合、それぞれの録画開始時刻を0秒として
時間binを揃えます。通常は、比較したい1録画の`*_imec0`を1回ずつ処理する
使い方を推奨します。

## 主なオプション

- `--bin-ms 100`: 時間bin幅。既定値は100 ms。
- `--min-spikes N`: スパイク数がN未満のgood clusterを除外。
- `--normalize sqrt_zscore`: 平方根変換後、クラスタごとにz-score。
- `--normalize zscore`: 生カウントを直接z-score。
- `--normalize none`: 正規化しない。
- `--max-neurons N`: 発火数上位Nクラスタだけを使用。0は全クラスタ。
- `--max-time-bins N`: 時間binを等間隔にN点へ間引く。0は全bin。
- `--sample-rate 30000`: `params.py`を使わずサンプルレートを指定。
- `--run-rastermap`: 行列作成後にRastermapも実行。
- `--n-pcs`、`--n-clusters`、`--locality`、`--time-lag-window`:
  Rastermapの設定。

全オプションは次で確認できます。

```powershell
python .\rastermap_good_clusters\run_good_clusters_rastermap.py --help
```

## 出力

`--out`で指定したフォルダに次を保存します。

- `spks_good.npy`: good clusterの発火数行列。
- `spks_good_normalized.npy`: Rastermap入力用の正規化行列。
- `time_seconds.npy`: 各時間binの時刻。
- `cluster_metadata.json`: 行番号、元データ、cluster ID、スパイク数。
- `summary.json`: データ数、good cluster数、bin数、処理条件。
- `good_clusters_unsorted.png`: Rastermap前のプレビュー。

`--run-rastermap`を付けた場合は、さらに次を保存します。

- `rastermap_result.npz`
- `rastermap_isort.npy`
- `cluster_metadata_rastermap_sorted.json`
- `good_clusters_rastermap.png`

## よくある問題

`No datasets with good clusters were found`

- 指定ルート以下に`imec*_ks4`があるか確認してください。
- Phyで保存した`cluster_group.tsv`に`good`があるか確認してください。
- フォルダ名が異なる場合は`--dataset-glob`を変更してください。

`rastermap is not installed`

```powershell
python -m pip install rastermap
```

メモリ不足の場合は、`--bin-ms`を大きくするか、
`--max-time-bins`、`--max-neurons`を指定してください。
