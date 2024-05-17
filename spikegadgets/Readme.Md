<root>
    <section title="Add trodesexport.exe to system path">
        <description>
            trodesexport.exeをシステムパスに追加する方法です。
            設定 → システムの詳細設定 → 環境変数 → システム環境編集 → 編集
        </description>
    </section>

    <section title="Spikegadgets: Concatinate REC Files">
        <command>
            trodesexport -kilosort -interp 0 -rec 20220930_145002.rec -rec 20220930_151129.rec -rec 20220930_153248.rec -rec 20220930_155359.rec
        </command>
        <description>
            上記のコマンドをCMDプロンプトで実行することで、複数のRECファイルを連結して.datファイルを作成できます。
        </description>
    </section>

    <section title="DIO Files">
        <command>
            trodesexport -dio -interp 0 -rec 20220930_145002.rec -rec 20220930_151129.rec -rec 20220930_153248.rec -rec 20220930_155359.rec
        </command>
        <description>
            DIOファイルを生成するコマンドです。
        </description>
    </section>

    <section title="Kilosort">
        <description>
            Gドライブの解析対象ディレクトリに移動し、getChanMapを実行してchanmap.matを作成します。???_channelmap_probe?.datファイルはchanmapフォルダを作成して移動させます（main_kilosort3がデータファイルと間違えるため）。
        </description>
    </section>

    <section title="Kilosort4">
        <description>
            メモリーがオーバーフローする場合は、Extra settingsでdeminxを32に設定します。CMDプロンプトで以下のようにOPENBLAS_NUM_THREADSを1に設定することも必要かもしれません。
        </description>
        <command>
            set OPENBLAS_NUM_THREADS=1
        </command>
    </section>

    <section title="Phy">
        <description>
            Anacondaでenvなしで起動し、以下のコマンドを実行します。
        </description>
        <command>
            phy template-gui params.py
        </command>
        <description>
            Ctrl+Alt+Nでノイズを指定し、Gでmergeします。Filter windowで「group!=’noise’」とすると見やすくなります。保存後、カレーション後（ks_postprocessingなど）、.phyフォルダをリネームします。
        </description>
    </section>

    <section title="Spikeinterface">
        <description>
            Anaconda上でsi_envというenvironmentにspikeinterfaceをインストールし、以下のように実行します。
        </description>
        <command>
            conda env si_env
        </command>
        <command>
            (si_env) C:\Users\stakahas.DESKTOP-JMVVEDC\Documents>jupyter notebook
        </command>
        <description>
            Vs codeでnp1_r1.ipynbをopenし、spike sortingはエラーになるので実行しません。stFolderと???.rawファイルが生成されます。
        </description>
    </section>

    <section title="Ecephys_spike_sorting Process">
        <description>
            Kilosort4を実行し、Kilosort4フォルダでphy template-gui params.pyを実行し、そのまま何もせずに保存します。その後、Vs codeでexportPhy.ipynbをopenし、base_folderに処理したいフォルダをセットします。mPhyフォルダが作成され、結果が入ります。
        </description>
        <steps>
            <step>
                <description>
                    phy template-gui params.pyでカレーションを行います。Alt+GでGood unitを設定して保存します。
                </description>
            </step>
            <step>
                <description>
                    再度、exportPhy.ipynbを実行します。remove_duplicated_spikesをms=1.4に変更し、waveforms-kilosort4をwaveforms-kilosort4-2に変更し、mPhyをmPhy2に変更します。
                </description>
            </step>
            <step>
                <description>
                    mPhy2に結果が入ります。mPhyはprePhyという名前に変更します（後処理でmPhyを自動検出してしまうため）。
                </description>
            </step>
            <step>
                <description>
                    最後に、mPhy2フォルダの中でphy template-gui params.pyを実行し、Alt+GですべてのユニットをGood unitに指定して終了します。
                </description>
            </step>
        </steps>
        <notes>
            <note>
                d-primeが~3.0以下はノイズの可能性あり。Similarity >0.85以上はmergeの可能性あり。Mono-synaptic or burstを考慮すべきです。
            </note>
        </notes>
    </section>

    <section title="DIO Data">
        <command>
            trodesexport -dio -rec 20220930_145002.rec
        </command>
        <description>
            MATLABコンソールで以下を実行します。
        </description>
        <matlab>
            <command>
                cd 20220930_145002.DIO
            </command>
            <command>
                fn='20230819_164611.dio_Controller_Din1.dat';
            </command>
            <command>
                a=readTrodesExtractedDataFile(fn);
            </command>
            <command>
                [b,c]=a.fields.data;
            </command>
        </matlab>
    </section>
</root>

