
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

Alt+Nでノイズを指定し、Gでmergeします。Filter windowで`group!=’noise’`とすると見やすくなります。

# Spikeinterface

Anaconda上で`si_env`というenvironmentにspikeinterfaceをインストールし、以下のように実行します。

### Commands

```sh
conda env si_env
(si_env) C:\Users\stakahas.DESKTOP-JMVVEDC\Documents>jupyter notebook
```

### Additional Info
jupyter notebookをVS code上で実行した方が便利かもしれません。その場合、si_env kernelを指定して実行する。
VS Codeで`np1_r1.ipynb`を実行し、`stFolder`と`???.raw`ファイルが生成されます。

# Ecephys_spike_sorting Process

Kilosort4を???.rawファイルを指定して実行し、Kilosort4フォルダで`phy template-gui params.py`を実行する。そのまま何もせずに保存します。その後、VS Codeで`exportPhy.ipynb`をopenし、`base_folder`に処理したいフォルダをセットします。`mPhy`フォルダが作成され、結果が入ります。

### Steps

1. mPhyフォルダ上で、`phy template-gui params.py`でcurationを行います。Alt+GでGood unitを設定して保存します。
2. `exportPhy2.ipynb`を実行します。exportPhy.ipynbとの違いは、`remove_duplicated_spikes`を`ms=1.4`に変更し、`waveforms-kilosort4`を`waveforms-kilosort4-2`に変更し、`mPhy`を`mPhy2`に変更したことです。
3. `mPhy2`に結果が入ります。`mPhy`は`prePhy`という名前に変更します（後処理で`mPhy`を自動検出してしまうため）。
4. 最後に、`mPhy2`フォルダの中で`phy template-gui params.py`を実行し、Alt+GですべてのユニットをGood unitに指定して終了します。

### Notes

- d-primeが~3.0以下はノイズの可能性あり。waveformを見てチェックしてください。不審な点がある場合は、`remove_duplicated_spikes`の`ms`パラメータを変更する。
- Similarity >0.85以上はmergeの可能性あり。Mono-synaptic or burstを考慮すべきです。

### MATLAB Commands
readSpks.mでデータを読み取る。

このスクリプト `readSpks.m` は、以下のディレクトリ構造を期待しています。

# ディレクトリ構造
## ルートディレクトリ
- `root_directory/`
  - `*.kilosort/`
    - `stFolder/`
      - `kilosort4/`
        - `mPhy*`
    - `timestamps.dat`
    - その他の.kilosortファイル
  - `*.DIO/`
    - `*dio?.dat`
  - `video/`
    - .matファイル（このフォルダには.matファイルが存在してはいけません）

## 詳細

1. **ルートディレクトリ**:
   - スクリプトの `directoryPath` 引数で指定されます。

2. **Kilosortディレクトリ**:
   - ルートディレクトリ内に1つ以上の`.kilosort`ディレクトリが存在する必要があります。
   - 各.kilosortディレクトリには `timestamps.dat` ファイルが含まれている必要があります。
   - `stFolder`ディレクトリが含まれ、その中に`kilosort4`ディレクトリと`mPhy*`ファイルが含まれている必要があります。

3. **DIOディレクトリ**:
   - ルートディレクトリ内に1つ以上の`.DIO`ディレクトリが存在する必要があります。
   - 各.DIOディレクトリには`*dio?.dat`ファイルが含まれている必要があります。

4. **videoディレクトリ**:
   - ルートディレクトリ内に`video`ディレクトリが存在し、.matファイルには,xydimsという変数名でxy座標が格納されている必要があります。



