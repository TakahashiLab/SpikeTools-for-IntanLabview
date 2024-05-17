
# Add trodesexport.exe to System Path

trodesexport.exeをシステムパスに追加する方法です。

- 設定 → システムの詳細設定 → 環境変数 → システム環境編集 → 編集

# Spikegadgets: Concatinate REC Files

### Command

```sh
trodesexport -kilosort -interp 0 -rec 20220930_145002.rec -rec 20220930_151129.rec -rec 20220930_153248.rec -rec 20220930_155359.rec
```

### Description

上記のコマンドをCMDプロンプトで実行することで、複数のRECファイルを連結して.datファイルを作成できます。

# DIO Files

### Command

```sh
trodesexport -dio -interp 0 -rec 20220930_145002.rec -rec 20220930_151129.rec -rec 20220930_153248.rec -rec 20220930_155359.rec
```

### Description

DIOファイルを生成するコマンドです。

# Kilosort

Gドライブの解析対象ディレクトリに移動し、`getChanMap`を実行して`chanmap.mat`を作成します。

`???_channelmap_probe?.dat`ファイルはchanmapフォルダを作成して移動させます（main_kilosort3がデータファイルと間違えるため）。

# Kilosort4

メモリーがオーバーフローする場合は、Extra settingsで`deminx`を32に設定します。

CMDプロンプトで以下のように`OPENBLAS_NUM_THREADS`を1に設定することも必要かもしれません。

### Command

```sh
set OPENBLAS_NUM_THREADS=1
```

# Phy

Anacondaでenvなしで起動し、以下のコマンドを実行します。

### Command

```sh
phy template-gui params.py
```

### Additional Info

Ctrl+Alt+Nでノイズを指定し、Gでmergeします。Filter windowで`group!=’noise’`とすると見やすくなります。保存後、カレーション後（ks_postprocessingなど）、`.phy`フォルダをリネームします。

# Spikeinterface

Anaconda上で`si_env`というenvironmentにspikeinterfaceをインストールし、以下のように実行します。

### Commands

```sh
conda env si_env
(si_env) C:\Users\stakahas.DESKTOP-JMVVEDC\Documents>jupyter notebook
```

### Additional Info

VS Codeで`np1_r1.ipynb`をopenし、spike sortingはエラーになるので実行しません。`stFolder`と`???.raw`ファイルが生成されます。

# Ecephys_spike_sorting Process

Kilosort4を実行し、Kilosort4フォルダで`phy template-gui params.py`を実行し、そのまま何もせずに保存します。その後、VS Codeで`exportPhy.ipynb`をopenし、`base_folder`に処理したいフォルダをセットします。`mPhy`フォルダが作成され、結果が入ります。

### Steps

1. `phy template-gui params.py`でカレーションを行います。Alt+GでGood unitを設定して保存します。
2. 再度、`exportPhy.ipynb`を実行します。`remove_duplicated_spikes`を`ms=1.4`に変更し、`waveforms-kilosort4`を`waveforms-kilosort4-2`に変更し、`mPhy`を`mPhy2`に変更します。
3. `mPhy2`に結果が入ります。`mPhy`は`prePhy`という名前に変更します（後処理で`mPhy`を自動検出してしまうため）。
4. 最後に、`mPhy2`フォルダの中で`phy template-gui params.py`を実行し、Alt+GですべてのユニットをGood unitに指定して終了します。

### Notes

- d-primeが~3.0以下はノイズの可能性あり。
- Similarity >0.85以上はmergeの可能性あり。Mono-synaptic or burstを考慮すべきです。

# DIO Data

### Command

```sh
trodesexport -dio -rec 20220930_145002.rec
```

### MATLAB Commands

```matlab
cd 20220930_145002.DIO
fn='20230819_164611.dio_Controller_Din1.dat';
a=readTrodesExtractedDataFile(fn);
[b,c]=a.fields.data;
```
