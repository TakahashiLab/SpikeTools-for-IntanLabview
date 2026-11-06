import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from prephy_feedback import FeedbackStore
from prephy_gui import build_command, inspect_datasets


class PrephyGUITest(unittest.TestCase):
    def test_commands_and_existing_model_protection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base = root / 'base model.joblib'
            base.touch()
            output = root / 'new model.joblib'
            command = build_command('train', [root], output, base, '2', '0.6')
            self.assertEqual(command[command.index('--input') + 1], str(root))
            self.assertEqual(command[command.index('--base-model') + 1], str(base))
            self.assertEqual(command[command.index('--threshold') + 1], '0.6')
            self.assertNotIn('--overwrite', command)
            with self.assertRaises(ValueError):
                build_command('train', [root], base)
            command = build_command('predict', [root], base)
            self.assertNotIn('--threshold', command)
            self.assertNotIn('--base-model', command)
            for weight in ['nan', 'inf', '-1']:
                with self.assertRaises(ValueError):
                    build_command('train', [root], output, weight=weight)

    def test_saved_feedback_and_changed_membership(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / 'imec1002_ks4'
            directory.mkdir()
            assignments = np.array([0, 0, 1, 1])
            np.save(directory / 'spike_clusters.npy', assignments)
            np.save(directory / 'spike_times.npy', np.arange(4))
            (directory / 'params.py').write_text('sample_rate = 30000\n')
            (directory / 'cluster_group.tsv').write_text('cluster_id\tgroup\n0\tgood\n1\tnoise\n')
            rows = inspect_datasets(tmp)
            self.assertEqual(rows[0][1:3], ('履歴なし', 0))
            store = FeedbackStore(directory)
            store.folder.mkdir()
            store.path.write_text(json.dumps({'reviewed': {}, 'save_number': 0}))
            store.commit(assignments, {0: 'good', 1: 'noise'}, complete=True)
            rows = inspect_datasets(tmp)
            self.assertEqual(rows[0][1:3], ('全確認完了', 2))
            before = store.path.read_bytes()
            inspect_datasets(tmp)
            self.assertEqual(store.path.read_bytes(), before)
            np.save(directory / 'spike_clusters.npy', np.array([0, 2, 1, 1]))
            rows = inspect_datasets(tmp)
            self.assertEqual(rows[0][1:3], ('保存後に変更あり', 0))


if __name__ == '__main__':
    unittest.main()
