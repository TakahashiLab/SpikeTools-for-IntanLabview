# Prephy GUI：別PCへの移植手順（Windows）

## 1. GitHubからコードを取得する

この手順書を含むブランチは `prephy-portable` です。GitHubへのアクセス権があるアカウントで取得してください。
PowerShellで、保存したいフォルダに移動して実行します。

```powershell
git clone --branch prephy-portable https://github.com/TakahashiLab/SpikeTools-for-IntanLabview.git
cd SpikeTools-for-IntanLabview\neuropixels
```

取得済みの場合は、未保存の変更を先にコミットするなどしてから、リポジトリ内で実行します。

```powershell
git fetch origin
git switch prephy-portable
git pull --ff-only origin prephy-portable
cd neuropixels
```

Gitはコードを移します。Python環境、学習済みモデル、解析データ、ユーザーの `.phy` 設定は別に用意します。
MATLABはPrephy GUIの起動・学習・予測には不要です。

## 2. Python環境を用意する

移植元で確認した構成は、64bit Python 3.11.3、NumPy 1.26.4、SciPy 1.10.1、
scikit-learn 1.2.2、joblib 1.2.0、threadpoolctl 2.2.0です。
既存モデルの読込には、特にscikit-learnのバージョンを揃えてください。

### GUIの学習・予測だけを使う場合

Tkinterを含む64bit Python 3.11を用意し、`neuropixels` 内で次を実行します。
以下の `py -3.11` は移植先に導入したPython 3.11を指します。

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-prephy.txt
.\.venv\Scripts\python.exe -c "import tkinter, numpy, scipy, sklearn, joblib; print(sklearn.__version__)"
```

最後に `1.2.2` と表示されることを確認します。Phy本体はこのrequirementsには含みません。

### 移植先に既存のPhy環境がある場合

既存環境のPythonを直接指定できます。環境に変更を加える前に、バージョンを確認してください。

```powershell
$env:PREPHY_PYTHON = 'C:\移植先のPython環境\python.exe'
& $env:PREPHY_PYTHON -c "import tkinter, numpy, scipy, sklearn, joblib; print(sklearn.__version__)"
.\start_prephy_gui.bat
```

`$env:PREPHY_PYTHON` の設定はそのPowerShell内だけ有効です。
デスクトップから起動する場合は、Windowsのユーザー環境変数に `PREPHY_PYTHON` を設定するか、
上の `.venv` を作成してください。起動BATは次の順でPythonを探します。

1. 環境変数 `PREPHY_PYTHON`（python.exeのフルパス）
2. `neuropixels\.venv\Scripts\python.exe`
3. `%USERPROFILE%\.conda\envs\phy2\python.exe`

移植元のPhy本体は `phy==2.1.0rc1`、`phylib==2.7.0` でした。
PhyのGUIには別途Qt等も必要です。上の最小環境はPrephy用で、Phy一式のインストーラーではありません。

## 3. 学習済みモデルと解析データをコピーする

使いたい `.joblib` と、存在すれば同名の `.joblib.json` を別PCにコピーします。
既定のモデル名は `neuropixels\models\good_model.joblib` です。
別の場所に置いてGUIの「参照…」から指定しても構いません。
環境変数 `PREPHY_MODEL` でも起動時のモデルを指定できます。
モデルはGitに同梱していません。

移植元GUIで設定されていた場所は `P:\seabird\models\good_model.joblib` でしたが、
移植作業時にはこちらからそのファイルの存在を確認できていません。元PCで実際に使っているモデルを確認してください。

予測だけなら、モデルと予測対象データが必要です。
再学習には教師データも必要で、`prephy_feedback` フォルダも解析結果と一緒にコピーしてください。
`.npy`、`params.py`、手動ラベルTSV、指標類を含む解析フォルダを保ち、
`params.py` 内のバイナリ参照先も移植先で読める場所に合わせます。

「元モデル」を使う再学習は、モデル内部に保存された元の教師データの絶対パスを再読込します。
ドライブ名やパスが変わると、その教師データを見つけられないことがあります。
元のパスを利用できるようにするか、「元モデル」を空欄にして、移植先で選択した確認済みデータから新規学習してください。
後者は元モデルの全教師データを自動継承する操作ではありません。

## 4. GUIを起動する

`start_prephy_gui.bat` をダブルクリックします。
検索フォルダは移植先のデータを選択してください。モデルの保存先も「参照…」で指定できます。
必要ならBATのショートカットをデスクトップに作成します。

起動だけを確認する場合はPowerShellから次を実行します。

```powershell
.\start_prephy_gui.bat --smoke-test
```

`GUI_SMOKE_OK` が表示されれば起動確認は成功です。学習・予測は実行しません。

## 5. Phyで確認履歴を記録する場合

GUIでの予測だけなら、この節のプラグイン導入は不要です。
Phyでの確認操作を記録して再学習に使う場合は、移植先のPhy環境にプラグインを設定します。
Phyを終了し、既存の `%USERPROFILE%\.phy` をバックアップしてから作業してください。

必要なプラグインを `%USERPROFILE%\.phy\plugins` にコピーします。

| ファイル | 用途 |
|---|---|
| `prephy_feedback.py` と `prephy_feedback_plugin.py` | 初回ラベル割り当て・確認履歴（再学習に必要） |
| `prephy_label_buffer_plugin.py` | 任意：good/noiseの一括操作 |
| `prephy_merge_candidates_plugin.py` | 任意：merge候補 |
| `similarity_left_columns.py` | 任意：Prephy列の表示順 |

`%USERPROFILE%\.phy\phy_config.py` の既存設定に、必要な名前を追加します。
既存のプラグインを残し、同じ名前は重複させないでください。
設定ファイルがまだない場合、次を最小例として作成できます。

```python
c = get_config()
c.TemplateGUI.plugins = [
    'PrephyFeedbackPlugin',
    'PrephyLabelBufferPlugin',
    'PrephyMergeCandidatesPlugin',
    'SimilarityLeftColumnsPlugin',
]
```

旧 `PhyGoodPolicyPlugin` は有効にしないでください。
merge候補は `.phy/good_policy_config.json` の設定を読むため、
元PCと同じ条件を使いたい場合はその設定もコピーします。

付属の `install_prephy_feedback.py`、`install_prephy_label_buffers.py`、
`install_prephy_merge_candidates.py` は、既存の `phy_config.py` に複数行の
`c.TemplateGUI.plugins = [...]` がある場合に利用できます。
それ以外の設定形式・新規環境では上記の手動設定を使ってください。
`install_prephy_columns.py --apply` は列表示のファイルだけを配置するため、
`SimilarityLeftColumnsPlugin` の設定への追加は別途必要です。

Phyを再起動し、確認操作を行ってCtrl+Sで保存した後、Phyを閉じてからPrephy GUIで検索します。
詳しい操作方法は [PREPHY_GOOD_README.md](PREPHY_GOOD_README.md) を参照してください。

## 6. 動作確認・よくある問題

最小環境での機能テスト（実データは変更しません）：

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s . -p "test_prephy*.py"
```

- Pythonが見つからない：BATの探索順と `PREPHY_PYTHON` を確認します。
- Tkinterがない：Tkinterを含むPythonを使います。`pip install tkinter` では導入しません。
- モデルが読めない：モデルとscikit-learn等のバージョンを確認します。
- 確認済みデータが出ない：Phyプラグイン設定と、保存済みの `prephy_feedback` を確認します。
- 予測ファイルが既にある：GUIは停止します。別の未予測データで確認してください。

この手順は移植元環境のコード・依存バージョンを基にしています。
移植先PCでのインストール・実データの学習／予測は別途確認が必要です。
