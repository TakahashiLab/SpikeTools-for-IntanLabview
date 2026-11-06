"""Meaningful end-to-end tests using disposable phy-format datasets."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

spec = importlib.util.spec_from_file_location('prephy_good', Path(__file__).with_name('prephy_good.py'))
app = importlib.util.module_from_spec(spec)
spec.loader.exec_module(app)


class WorkflowTests(unittest.TestCase):
    def make_data(self, root, recording, shank):
        directory = root / ('catgt_' + recording + '_g0') / (recording + '_g0_imec0') / ('imec' + str(shank) + '_ks4')
        directory.mkdir(parents=True)
        (directory / 'params.py').write_text("sample_rate = 1000\ndat_path = '../missing.ap.bin'\ndtype = 'int16'\nn_channels_dat = 4\n", encoding='utf-8')
        np.save(directory / 'spike_clusters.npy', np.repeat(np.arange(20), 12))
        np.save(directory / 'spike_times.npy', np.arange(240) * 10)
        (directory / 'cluster_group.tsv').write_text('cluster_id\tgroup\n' + ''.join(f'{i}\t{"good" if i < 10 else "noise"}\n' for i in range(20)) + '999\tgood\n', encoding='utf-8')
        # Opposite automatic labels must not replace manual curation.
        (directory / 'cluster_KSLabel.tsv').write_text('cluster_id\tKSLabel\n' + ''.join(f'{i}\t{"mua" if i < 10 else "good"}\n' for i in range(20)), encoding='utf-8')
        (directory / 'metrics.csv').write_text('cluster_id,amplitude,snr,isi_viol\n' + ''.join(f'{i},{100 if i < 10 else 5},{10 if i < 10 else 0.1},{0 if i < 10 else 1}\n' for i in range(20)), encoding='utf-8')
        return directory

    def test_train_predict_preserves_curation_and_separates_recordings(self):
        from argparse import Namespace
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'curated'
            for run in ('A', 'B', 'C', 'D'):
                for shank in (1000, 1001):
                    self.make_data(source, run, shank)
            excluded = source / 'catgt_A_g0/A_g0_imec0/imec1000_ks4_orig'
            excluded.mkdir()
            (excluded / 'params.py').write_text('')
            np.save(excluded / 'spike_clusters.npy', [0])
            model = root / 'model.joblib'
            app.train(Namespace(input=[source], model=model, overwrite=False,
                                validation_fraction=0.25, threshold=0.5, jobs=1))
            metadata = json.loads(model.with_suffix('.joblib.json').read_text())
            self.assertEqual(metadata['training_clusters'], 160)
            self.assertEqual(metadata['training_good'], 80)
            validation = metadata['validation']
            self.assertEqual(validation['status'], 'evaluated')
            self.assertFalse(set(validation['training_recordings']) & set(validation['heldout_recordings']))
            self.assertGreater(validation['precision'], 0.9)
            target = self.make_data(root / 'new', 'E', 1000)
            protected = [target / 'cluster_group.tsv', target / 'spike_clusters.npy', target / 'cluster_KSLabel.tsv']
            before = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in protected}
            options = Namespace(input=[target], model=model, overwrite=False, threshold=None,
                                min_feature_fraction=0.5, report=root / 'report.json')
            app.predict(options)
            rows = app.read_rows(target / 'prephy_good_predictions.tsv')
            self.assertEqual(len(rows), 20)  # Ignore deleted/inactive label 999.
            self.assertEqual(rows[0]['prephy_good_prediction'], 'good')
            self.assertEqual(rows[-1]['prephy_good_prediction'], 'not_good')
            for path in protected:
                self.assertEqual(before[path], hashlib.sha256(path.read_bytes()).hexdigest())
            with self.assertRaises(ValueError):
                app.predict(options)
            # Only a saved explicit correction is added when retraining from a
            # base model, rather than all automatically predicted labels.
            from prephy_feedback import FeedbackStore
            store = FeedbackStore(target)
            assignments = np.load(target / 'spike_clusters.npy')
            store.initialize(assignments)
            store.initialized()
            labels = app.manual_labels(target)
            labels[0] = 'noise'
            app.write_tsv(target / 'cluster_group.tsv', ['cluster_id', 'group'], sorted(labels.items()))
            store.commit(assignments, labels, reviewed_ids=[0])
            updated_model = root / 'updated.joblib'
            app.train(Namespace(input=[target], base_model=model, model=updated_model,
                                overwrite=False, validation_fraction=0.25,
                                threshold=0.5, jobs=1, feedback_weight=2.0))
            updated = json.loads(updated_model.with_suffix('.joblib.json').read_text())
            self.assertEqual(updated['training_clusters'], 161)
            feedback_sources = [s for s in updated['datasets'] if s['label_source'] == 'saved_explicit_feedback']
            self.assertEqual(len(feedback_sources), 1)
            self.assertEqual(feedback_sources[0]['labels'], {'noise': 1})
            (target / 'metrics.csv').unlink()
            options.overwrite = True
            app.predict(options)
            self.assertTrue(all(row['prephy_good_prediction'] == 'unknown' for row in app.read_rows(target / 'prephy_good_predictions.tsv')))

    def test_no_kslabel_fallback_and_no_execution_of_params(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = self.make_data(Path(temporary), 'A', 1000)
            (directory / 'cluster_group.tsv').unlink()
            self.assertEqual(app.manual_labels(directory), {})
            marker = directory / 'executed.txt'
            with (directory / 'params.py').open('a') as handle:
                handle.write(f"\nopen({str(marker)!r}, 'w').write('bad')\n")
            ids, _, _, _ = app.extract(directory)
            self.assertEqual(len(ids), 20)
            self.assertFalse(marker.exists())


if __name__ == '__main__':
    unittest.main()
