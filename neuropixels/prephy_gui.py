"""Desktop interface for saved Prephy feedback, retraining and prediction."""
import json
import os
from pathlib import Path
import queue
import subprocess
import sys
import threading
from datetime import datetime
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

HERE = Path(__file__).resolve().parent
MODEL = Path(os.environ.get('PREPHY_MODEL', str(HERE / 'models' / 'good_model.joblib'))).expanduser()


def inspect_datasets(root):
    """Read feedback status; never initialize or modify a dataset."""
    from prephy_good import discover
    from prephy_feedback import training_feedback
    rows = []
    for directory in discover([Path(root)]):
        state_path = directory / 'prephy_feedback' / 'state.json'
        try:
            if not state_path.is_file():
                rows.append((str(directory), '履歴なし', 0, ''))
                continue
            state = json.loads(state_path.read_text(encoding='utf-8'))
            usable = training_feedback(directory)
            recorded = len(state.get('reviewed', {}))
            valid = len(usable or {})
            if not state.get('saved_membership_sha256'):
                status = '未保存'
            elif valid != recorded:
                status = '保存後に変更あり'
            elif not valid:
                status = '未確認'
            else:
                status = '全確認完了' if state.get('complete') else '一部確認済み'
            rows.append((str(directory), status, valid, state.get('updated', '')))
        except Exception as error:
            rows.append((str(directory), '読込エラー', 0, str(error)))
    return rows


def build_command(mode, inputs, model, base='', weight='2', threshold=''):
    model = Path(model).expanduser()
    if model.suffix.lower() != '.joblib':
        raise ValueError('モデルは .joblib ファイルを指定してください。')
    if not inputs:
        raise ValueError('対象フォルダを指定してください。')
    for item in inputs:
        if not Path(item).is_dir():
            raise ValueError('フォルダがありません: ' + str(item))
    command = [sys.executable, '-u', str(HERE / 'prephy_good.py'), mode]
    for item in inputs:
        command.extend(['--input', str(item)])
    command.extend(['--model', str(model)])
    if mode == 'train':
        if model.exists():
            raise ValueError('出力モデルが既にあります。別の名前を指定してください。')
        if not model.parent.is_dir():
            raise ValueError('出力先フォルダがありません。')
        if base:
            if not Path(base).is_file():
                raise ValueError('元モデルがありません。')
            command.extend(['--base-model', str(base)])
        import math
        if not math.isfinite(float(weight)) or float(weight) <= 0:
            raise ValueError('教師データの重みは正の数を指定してください。')
        command.extend(['--feedback-weight', str(weight)])
    elif mode == 'predict':
        if not model.is_file():
            raise ValueError('予測モデルがありません。')
    else:
        raise ValueError('Unknown mode')
    if threshold.strip():
        if not 0 < float(threshold) < 1:
            raise ValueError('閾値は0より大きく1より小さい値にしてください。')
        command.extend(['--threshold', threshold])
    return command


class PrephyGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('Prephy — 確認済みデータ・再学習・予測')
        self.geometry('1120x800')
        self.minsize(900, 680)
        self.events = queue.Queue()
        self.busy = False
        self.rows = {}
        self.process = None
        self.status = tk.StringVar(value='フォルダを検索し、学習に使う確認済みデータを選択してください。')
        self.protocol('WM_DELETE_WINDOW', self.close)
        style = ttk.Style(self)
        style.theme_use('clam')
        style.configure('Treeview', rowheight=28)
        style.configure('TButton', padding=5)
        outer = ttk.Frame(self, padding=14)
        outer.pack(fill='both', expand=True)
        ttk.Label(outer, text='Prephy', font=('Yu Gothic UI', 20, 'bold')).pack(anchor='w')
        ttk.Label(outer, text='Phyで確認・保存 → データを選択 → 再学習 → 次のデータを予測').pack(anchor='w', pady=(0, 10))
        tabs = ttk.Notebook(outer)
        tabs.pack(fill='both', expand=True)
        training = ttk.Frame(tabs, padding=12)
        prediction = ttk.Frame(tabs, padding=12)
        tabs.add(training, text='確認済みデータ / 再学習')
        tabs.add(prediction, text='予測')
        self.search = self.path_row(training, '検索フォルダ', '', 'dir')
        controls = ttk.Frame(training)
        controls.pack(fill='x', pady=5)
        self.scan_button = ttk.Button(controls, text='保存履歴を検索 / 更新', command=self.scan)
        self.scan_button.pack(side='left')
        ttk.Button(controls, text='全確認完了をまとめて選択', command=self.select_complete).pack(side='left', padx=6)
        ttk.Label(training, text='Ctrl / Shiftで複数選択。「確認済み数」は現在の保存データと照合した件数です。').pack(anchor='w')
        table_frame = ttk.Frame(training)
        table_frame.pack(fill='both', expand=True, pady=6)
        self.table = ttk.Treeview(table_frame, columns=('status', 'count', 'date', 'path'), show='headings', selectmode='extended', height=7)
        for col, title, width in [('status', '確認状態', 135), ('count', '確認済み数', 85), ('date', '最終保存（UTC）', 185), ('path', '解析フォルダ', 600)]:
            self.table.heading(col, text=title)
            self.table.column(col, width=width, minwidth=70, stretch=col == 'path')
        self.table.grid(row=0, column=0, sticky='nsew')
        yscroll = ttk.Scrollbar(table_frame, orient='vertical', command=self.table.yview)
        yscroll.grid(row=0, column=1, sticky='ns')
        xscroll = ttk.Scrollbar(table_frame, orient='horizontal', command=self.table.xview)
        xscroll.grid(row=1, column=0, sticky='ew')
        self.table.configure(yscrollcommand=yscroll.set, xscrollcommand=xscroll.set)
        table_frame.rowconfigure(0, weight=1)
        table_frame.columnconfigure(0, weight=1)
        self.base = self.path_row(training, '元モデル', str(MODEL), 'open')
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        self.output = self.path_row(training, '新モデルの保存先', str(MODEL.with_name('good_model_' + stamp + '.joblib')), 'save')
        options = ttk.Frame(training)
        options.pack(fill='x', pady=6)
        ttk.Label(options, text='手動確認データの重み').pack(side='left')
        self.weight = tk.StringVar(value='2')
        ttk.Entry(options, textvariable=self.weight, width=6).pack(side='left', padx=8)
        ttk.Label(options, text='good閾値').pack(side='left')
        self.train_threshold = tk.StringVar(value='0.5')
        ttk.Entry(options, textvariable=self.train_threshold, width=6).pack(side='left', padx=8)
        self.train_button = ttk.Button(training, text='選択した確認済みデータで再学習', command=self.train)
        self.train_button.pack(anchor='e')
        ttk.Label(training, text='Phyを保存して終了してから実行してください。元モデルに記録された教師データも再読込します。').pack(anchor='w', pady=5)
        self.predict_input = self.path_row(prediction, '予測するフォルダ', '', 'dir')
        self.predict_model = self.path_row(prediction, '予測モデル', str(MODEL), 'open')
        opts = ttk.Frame(prediction)
        opts.pack(fill='x', pady=12)
        ttk.Label(opts, text='good閾値（空欄ならモデルの値）').pack(side='left')
        self.predict_threshold = tk.StringVar()
        ttk.Entry(opts, textvariable=self.predict_threshold, width=8).pack(side='left', padx=8)
        ttk.Label(prediction, text='予測完了後にPhyで開いてください。保存済みの手動ラベルは保護されます。\n既存の予測ファイルがある場合は停止します。').pack(anchor='w', pady=8)
        self.predict_button = ttk.Button(prediction, text='予測ファイルを作成', command=self.predict)
        self.predict_button.pack(anchor='e')
        ttk.Label(outer, textvariable=self.status, wraplength=1050).pack(anchor='w', pady=(10, 5))
        log_frame = ttk.LabelFrame(outer, text='実行ログ', padding=5)
        log_frame.pack(fill='both')
        self.log = tk.Text(log_frame, height=8, wrap='word', state='disabled', font=('Consolas', 10))
        self.log.pack(fill='both', expand=True)
        self.after(100, self.poll)

    def path_row(self, parent, title, default, kind):
        row = ttk.Frame(parent)
        row.pack(fill='x', pady=4)
        ttk.Label(row, text=title, width=19).pack(side='left')
        value = tk.StringVar(value=default)
        ttk.Entry(row, textvariable=value).pack(side='left', fill='x', expand=True)
        def browse():
            if kind == 'dir':
                result = filedialog.askdirectory(parent=self)
            elif kind == 'save':
                result = filedialog.asksaveasfilename(parent=self, defaultextension='.joblib', filetypes=[('Prephy model', '*.joblib')])
            else:
                result = filedialog.askopenfilename(parent=self, filetypes=[('Prephy model', '*.joblib')])
            if result:
                value.set(result)
        ttk.Button(row, text='参照…', command=browse).pack(side='left', padx=(6, 0))
        return value

    def set_busy(self, value):
        self.busy = value
        for button in (self.scan_button, self.train_button, self.predict_button):
            button.configure(state='disabled' if value else 'normal')

    def scan(self):
        if self.busy:
            return
        root = self.search.get().strip()
        if not root or not Path(root).is_dir():
            messagebox.showerror('検索できません', '存在する検索フォルダを指定してください。', parent=self)
            return
        self.set_busy(True)
        self.status.set('保存履歴を検索・照合しています…')
        def work():
            try:
                self.events.put(('scan', inspect_datasets(root)))
            except Exception as error:
                self.events.put(('error', str(error)))
        threading.Thread(target=work, daemon=True).start()

    def select_complete(self):
        self.table.selection_set([key for key, row in self.rows.items() if row[1] == '全確認完了' and row[2] > 0])

    def train(self):
        if self.busy:
            return
        selected = [self.rows[key] for key in self.table.selection()]
        if not selected or any(row[2] <= 0 for row in selected):
            messagebox.showerror('再学習できません', '確認済み数が1以上のデータだけを選択してください。', parent=self)
            return
        self.start('train', [r[0] for r in selected], self.output.get().strip(),
                   base=self.base.get().strip(), weight=self.weight.get(), threshold=self.train_threshold.get())

    def predict(self):
        if not self.busy:
            self.start('predict', [self.predict_input.get().strip()], self.predict_model.get().strip(), threshold=self.predict_threshold.get())

    def start(self, mode, inputs, model, **options):
        try:
            if any(not p for p in inputs):
                raise ValueError('対象フォルダを指定してください。')
            command = build_command(mode, inputs, model, **options)
        except (ValueError, OSError) as error:
            messagebox.showerror('入力を確認してください', str(error), parent=self)
            return
        self.set_busy(True)
        self.status.set('再学習中…' if mode == 'train' else '予測中…')
        self.append_log('\n' + subprocess.list2cmdline(command) + '\n')
        def work():
            try:
                env = dict(os.environ, PYTHONIOENCODING='utf-8')
                self.process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, encoding='utf-8', errors='replace', env=env,
                    creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
                for line in self.process.stdout:
                    self.events.put(('log', line))
                code = self.process.wait()
                self.events.put(('done', (mode, code, model)))
            except Exception as error:
                self.events.put(('error', str(error)))
        threading.Thread(target=work, daemon=True).start()

    def append_log(self, text):
        self.log.configure(state='normal')
        self.log.insert('end', text)
        self.log.see('end')
        self.log.configure(state='disabled')

    def poll(self):
        try:
            while True:
                kind, value = self.events.get_nowait()
                if kind == 'log':
                    self.append_log(value)
                elif kind == 'scan':
                    self.table.delete(*self.table.get_children())
                    self.rows.clear()
                    for row in value:
                        key = self.table.insert('', 'end', values=(row[1], row[2], row[3], row[0]))
                        self.rows[key] = row
                    self.set_busy(False)
                    self.status.set('%dフォルダを確認。学習に使う行を選択してください。' % len(value))
                elif kind == 'done':
                    mode, code, model = value
                    self.set_busy(False)
                    self.process = None
                    if code == 0:
                        self.status.set('完了: ' + model if mode == 'train' else '予測完了。対象データをPhyで開いてください。')
                        if mode == 'train':
                            self.predict_model.set(model)
                    else:
                        self.status.set('実行失敗（終了コード%d）。ログを確認してください。' % code)
                elif kind == 'error':
                    self.set_busy(False)
                    self.status.set('処理に失敗しました。ログを確認してください。')
                    self.append_log(value + '\n')
        except queue.Empty:
            pass
        self.after(100, self.poll)

    def close(self):
        if self.busy:
            messagebox.showinfo('処理中', '処理が終了してから閉じてください。', parent=self)
        else:
            self.destroy()


if __name__ == '__main__':
    app = PrephyGUI()
    if '--smoke-test' in sys.argv:
        app.update()
        assert app.winfo_width() >= 900
        assert len(app.table['columns']) == 4
        app.destroy()
        print('GUI_SMOKE_OK')
    else:
        app.mainloop()
