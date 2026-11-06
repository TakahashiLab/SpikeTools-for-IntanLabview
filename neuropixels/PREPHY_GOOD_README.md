# phy起動前のgood判定

別PCへの導入は [PREPHY_MIGRATION.md](PREPHY_MIGRATION.md) を参照してください。

## GUIで確認済みデータの検索・再学習・予測

`start_prephy_gui.bat`をダブルクリックして起動します。GUIはPython標準のTkinterを使います。起動時は環境変数`PREPHY_PYTHON`、同じフォルダの`.venv`、ユーザーフォルダの`.conda/envs/phy2`の順でPythonを探します。モデルの初期値は`models/good_model.joblib`で、環境変数`PREPHY_MODEL`でも指定できます。

1. Phyで確認操作（選択確認または全確認完了）と保存を済ませ、Phyを終了します。
2. GUIの「検索フォルダ」に解析フォルダ、catgt、または個体のフォルダを指定し、「保存履歴を検索 / 更新」を押します。
3. 一覧で学習に使う行を選択します。Ctrl / Shiftで複数選択できます。「全確認完了をまとめて選択」も使えます。一部確認済みのデータも、確認済み数が1以上なら選択できます。
4. 元モデルと新モデルの保存先を指定して「選択した確認済みデータで再学習」を押します。元モデルを空欄にすると選択データだけで新規学習します。同名モデルの上書きは行いません。
5. 「予測」タブで予測先とモデルを指定して実行します。再学習が成功すると予測モデル欄に新モデルが入ります。予測後にPhyを起動してください。

確認済み数は保存された所属配列・ラベルと照合します。履歴なし・未確認のデータはGUIからの再学習に選べません。履歴のない過去データを教師にする場合は、Phyで明示的に確認・保存するか、下記のCLIを使います。保存後に変更があるデータはPhyで再度確認してください。読み取りエラーは一覧に表示します。

GUIは既存の`prephy_good.py`を実行し、ログと成功・失敗を表示します。予測は既存の手動ラベルを変更せず、既存の予測ファイルがある場合は停止します。学習や予測はボタンを押すまで実行しません。処理中はウィンドウを閉じられません。

`prephy_good.py` は手動修正済みのphyデータでRandom Forestを学習し、指定フォルダ以下のデータにgoodスコアを付けます。GUIは起動しません。`cluster_group.tsv`、スパイク所属、既存の手動修正には書き込みません。

バージョン1.1では、別の軽量プラグイン`PrephyFeedbackPlugin`が、phy初回起動時のラベル割り当てと保存済み手動修正の記録を担当します。旧`PhyGoodPolicyPlugin`は無効のままです。

## 実行例

PowerShellから、動作確認済みの`phy2`のPythonとスクリプトをフルパスで指定します。現在のフォルダやcondaの有効化状態に依存しません。依存関係はnumpy、scikit-learn、joblibです。

```powershell
$phyPython = "C:\Users\stakahas.DESKTOP-JMVVEDC\.conda\envs\phy2\python.exe"
$prephyScript = "J:\OneDrive - 同志社大学\Documents\MATLAB\SpikeTools-for-IntanLabview\neuropixels\prephy_good.py"
& $phyPython $prephyScript train --input "T:\seabird\2025\original\subjectT\20250906\catgt" --model "P:\seabird\models\good_model.joblib"
& $phyPython $prephyScript predict --input "P:\seabird\2026\subjectB1" --model "P:\seabird\models\good_model.joblib" --report "P:\seabird\models\B1_predictions.json"
```

上の例は学習と判定を別々に実行するものです。初回は`train`を実行し、その終了後、判定対象の解析完了を確認して`predict`を実行してください。`--input`は複数回指定でき、個別の`imec*_ks4`でも、その上位ディレクトリでも使えます。`_orig`、`_failed`、`.phy`の中は探索しません。処理中のデータには実行せず、後処理・手動保存の完了したフォルダを指定してください。

`C:\ProgramData\anaconda3\lib\site-packages\scipy`のDLLエラーが出る場合、base環境のPythonを使っています。上記の`$phyPython`で起動してください。`& $phyPython $prephyScript --version`が`1.1.0`を表示すれば、必要なライブラリを読み込めています。

## 教師ラベルと特徴量

- 教師ラベルは**`cluster_group.tsv`だけ**から読みます。goodを正例、mua/noiseを負例とし、unsorted、未判定、現在スパイクが所属しない削除済みクラスタは除外します。`cluster_KSLabel.tsv`での上書き・代用はしません。
- 教師データには自分で手動確認済みのフォルダを指定してください。ファイルの存在だけでは手動確認済みかは判別できません。過去にプラグインで自動付与したラベルも含まれ得ます。
- 特徴量はスパイク数、発火率、振幅、ContamPct、ISI違反、テンプレート最大類似度、波形幅、PT比、波形勾配、SNR、広がり、伝播速度です。`cluster_info.tsv`だけでなく、phy起動前に存在する`metrics.csv`や各種TSV/波形指標も読みます。
- スパイク数・発火率は現在の`spike_clusters.npy`から計算します。記録長はメタデータ→バイナリサイズ→最終スパイクの順で取得し、採用元を記録します。`params.py`は実行せず、リテラルの設定値だけを読みます。
- 電極番号・深さ・クラスタIDは学習特徴に使いません。異なるシャンク・個体の位置を品質として学習することを避けます。
- プラグインを参考にした独立モデルです。追加の波形特徴や特徴欠損を使用するため、既存プラグインのスコア・不応期ペナルティと数値的に同一ではありません。

## 学習と評価

記録の25%を評価用に分け、同じ記録の全シャンクを同じ側に置いてprecision、recall、混同行列、ROC-AUC、average precisionを計算します。欠損補完値も学習側だけで求めます。最後に全教師データでモデルを学習し直します。1記録しかない場合などは、独立評価不能として報告します。

モデルの隣に`.joblib.json`を保存します。教師データのパス・ラベル数・入力ファイル更新日時、特徴量、評価結果、バージョンが含まれます。同じ個体の記録が学習・評価の両方に含まれる可能性があるため、新しい個体や2026年データでの精度を保証する評価ではありません。

判定閾値は既定で0.5です。プラグイン設定の0.1を引き継ぎません。`train --threshold 0.8`または`predict --threshold 0.8`で変更できます。スコアはモデルの推定値で、精度が実測された確率を保証するものではありません。

## 出力

各データセットに以下の新規ファイルを保存します。

| ファイル | 内容 |
|---|---|
| `prephy_good_predictions.tsv` | スコア、good/not_good/unknown、特徴量充足率、元のgroup |
| `cluster_prephy_good_score.tsv` | phyで読める数値メタデータ列 `prephy_good_score` |
| `cluster_prephy_good_prediction.tsv` | phyで読める分類メタデータ列 `prephy_good_prediction` |
| `prephy_good_prediction.json` | モデルのハッシュ、閾値、入力ファイル情報、件数 |

学習に使った特徴の半数未満しか揃わないクラスタは`unknown`とし、phy用スコア列を空欄にします。`--min-feature-fraction`で変更できます。`not_good`はmua/noiseの区別を意味しません。既存出力があれば停止します。`--overwrite`は予測ファイル、またはtrain時のモデルを置き換えるオプションです。手動ラベルを変更するオプションはありません。

列表示用の`SimilarityLeftColumnsPlugin`はprephy列をID直後に表示するよう設定済みです。旧`PhyGoodPolicyPlugin`は無効にし、初回割り当て・記録には以下の新プラグインだけを使います。

## phyでの初回割り当てと手動修正の記録

`predict`が出力したデータを、再起動したphyで開きます。`PrephyFeedbackPlugin`は初回だけ、予測`good`を`good`、`not_good`と`unknown`を`noise`に割り当てます。**Ctrl+Sで保存してください。** 手動groupラベルが既にあるデータには初回割り当てを行いません。クラスタ構成が予測時と異なる場合も割り当てを中止します。

初回割り当てはクラスタID・所属の情報を照合し、各データセットの`prephy_feedback/state.json`に記録します。以後は再起動やpredictの再実行で手動ラベルを自動上書きしません。バージョン1.0の予測ファイルも、予測時のファイル更新日時・サイズが一致すれば使えます。不一致の場合はpredictを再実行してください。

| phyでの操作 | 記録・学習への扱い |
|---|---|
| good/noise/muaを手動変更してCtrl+S | そのクラスタを明示的な確認済み教師データとして記録 |
| `Alt+Shift+V` または `:prephyreview` | 選択クラスタは予測どおりで正しいと確認し、保存 |
| `Alt+Shift+C` または `:prephycomplete` | 全クラスタを確認済みとして保存（全てgood/noise/muaへの割り当てが必要） |
| 自動割り当てを保存しただけ | 未確認のまま。教師データには追加しない |
| merge/splitして保存 | 追加・削除ID、所属配列、特徴量を記録。新クラスタは明示的確認まで教師データにしない |

操作は`Tools → Prephy feedback`からも選べます。完了操作は実際に全クラスタを確認した後に使ってください。未保存の操作は学習に入りません。

### good/noiseの手動バッファ

`PrephyLabelBufferPlugin`で次の6機能を利用できます。旧判定プラグインの学習・候補提示は起動しません。merge候補だけの独立プラグインは下記を参照してください。予測ファイルのないデータでも利用できます。

| 操作 | コマンド | ショートカット |
|---|---|---|
| goodバッファに追加 | `:addgood` | `Alt+Shift+J` |
| goodバッファから除外 | `:unbuffergood` | `Alt+Shift+K` |
| goodバッファを一括適用 | `:applygood` | `Alt+Shift+L` |
| noiseバッファに追加 | `:addnoise` | `Alt+Shift+A` |
| noiseバッファから除外 | `:unbuffernoise` | `Alt+Shift+U` |
| noiseバッファを一括適用 | `:applynoise` | `Alt+Shift+B` |

`Tools → Label buffers`からも選べます。バッファ追加・除外ではラベルは変わりません。反対側に追加すると元のバッファから外れます。適用時にまとめて変更し、Ctrl+Sで保存します。フィードバック対象のデータでは、その保存時に手動修正として記録します。バッファはメモリ上だけで、phyを閉じると消えます。既に目的のラベルになっているクラスタは適用時にスキップします。そのようなクラスタを学習用に確認済みとするには、別途`Alt+Shift+V`を使用してください。

記録先はデータセット直下の`prephy_feedback/`です。初回のラベル・予測・モデル情報、保存ごとの操作と最終group、特徴量（cluster_info.tsv）、クラスタ所属のスナップショットを保持します。merge/split後は所属の変わったクラスタの確認状態を無効にし、変更のないクラスタの確認状態は保持します。このフォルダは学習履歴なので削除しないでください。

## merge候補だけを表示するプラグイン

`PrephyMergeCandidatesPlugin`は旧プラグインのmerge候補操作だけを独立させたものです。学習、モデル読込、自動ラベル割り当ては行いません。既存のフィードバックとラベルバッファは引き続き使用できます。

| 操作 | コマンド | ショートカット |
|---|---|---|
| 次の候補ペアを選択 | `:nextmerge` | `Alt+Shift+M` |
| 選択クラスタの候補を最大20個選択 | `:mergesuggest 20` | — |
| 現在選択中の2クラスタをmergeし次の候補へ | `:domerge` | `Alt+Shift+Enter` |
| 現在のペアを却下し次の候補へ | `:rejectmerge` | `Ctrl+Shift+R` |

`Tools → Merge candidates`からも操作できます。候補表示だけではmergeしません。merge後はCtrl+Sで保存してください。学習用の確認・再学習は前述のフィードバック操作で行います。却下したペアはその起動中だけ除外し、再起動するとリセットします。

候補はnoise以外から、現在のクラスタ所属に対応するテンプレート類似度と交差不応期で選びます。`~/.phy/good_policy_config.json`の`merge_similarity_threshold`（現在0.81）と`xcorr_refractory_ms`（現在1.2 ms）だけを読み、±不応期内に交差スパイクが1個でもあれば除外します。merge/split後の新IDにも対応します。旧プラグインの学習スコアを使うfallback判定は含みません。候補がない場合は条件に合うペアがないと表示します。

導入は`phy2`のPythonで`install_prephy_merge_candidates.py`を実行し、Phyを再起動します。設定と既存プラグインは変更前にバックアップします。旧`PhyGoodPolicyPlugin`は無効のままにしてください。

## 手動修正を反映した再学習

phyで保存して閉じた後、元モデルの教師データに、確認済みの手動修正を加えてモデルを作り直します。

```powershell
$phyPython = "C:\Users\stakahas.DESKTOP-JMVVEDC\.conda\envs\phy2\python.exe"
$prephyScript = "J:\OneDrive - 同志社大学\Documents\MATLAB\SpikeTools-for-IntanLabview\neuropixels\prephy_good.py"
& $phyPython $prephyScript train --base-model "P:\seabird\models\good_model.joblib" --input "T:\seabird\2025\original\subjectX\20250920\catgt" --model "P:\seabird\models\good_model_v2.joblib"
& $phyPython $prephyScript predict --input "次に判定するフォルダ" --model "P:\seabird\models\good_model_v2.joblib"
```

`--base-model`を指定すると、元モデルに記録された教師データの場所を再読込し、`--input`と合わせて再学習します。元モデルは保持します。フィードバック記録のあるデータは、保存済み・確認済みクラスタだけを採用します。自動割り当て全体を教師にすることはありません。手動フィードバックの重みは既定で2倍です（`--feedback-weight 1`で通常と同じ重み）。モデルの更新はこの再学習を実行した時点で反映され、phyの保存だけではモデルは変わりません。

フィードバック導入前の手動修正済みフォルダについては、従来どおりユーザー指定の`cluster_group.tsv`を教師とします。既にラベルを変更していたデータを後から開いて記録を開始した場合、過去の確認状態は自動推測しないため、確認済みクラスタの選択確認または全確認完了の操作を行ってください。

モデルは同じscikit-learnバージョンで読み込んでください。joblibモデルは自分で作成した信頼できるファイルだけを指定してください。

## 検証

```powershell
python test_prephy_good.py
```

手動ラベル優先、記録単位の学習・評価分離、予測TSV、既存ラベル・スパイク配列の保護、特徴不足での保留、params.pyを実行しないことを検証します。
