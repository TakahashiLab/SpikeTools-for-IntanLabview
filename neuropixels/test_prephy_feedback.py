from pathlib import Path
import tempfile
import unittest
import json

import numpy as np
from prephy_feedback import FeedbackStore, array_hash, training_feedback


class FeedbackTests(unittest.TestCase):
    def setup_dataset(self, path):
        assignments = np.repeat([0, 1, 2], 4)
        np.save(path / 'spike_clusters.npy', assignments)
        (path / 'cluster_group.tsv').write_text('cluster_id\tKSLabel\n0\tgood\n1\tmua\n2\tgood\n')
        (path / 'prephy_good_predictions.tsv').write_text('cluster_id\tprephy_good_prediction\n0\tgood\n1\tnot_good\n2\tunknown\n')
        (path / 'prephy_good_prediction.json').write_text(json.dumps({'membership_sha256': array_hash(assignments)}))
        return assignments

    def save_groups(self, path, groups):
        (path / 'cluster_group.tsv').write_text('cluster_id\tgroup\n' + ''.join(f'{cid}\t{label}\n' for cid, label in groups.items()))

    def test_auto_labels_are_not_teachers_and_manual_save_is(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            assignments = self.setup_dataset(path)
            store = FeedbackStore(path)
            self.assertEqual(store.initialize(assignments), {0: 'good', 1: 'noise', 2: 'noise'})
            store.initialized()
            self.assertEqual(store.initialize(assignments), {0: 'good', 1: 'noise', 2: 'noise'})
            altered = assignments.copy()
            altered[0] = 1
            with self.assertRaises(ValueError):
                store.initialize(altered)
            groups = {0: 'good', 1: 'noise', 2: 'noise'}
            self.save_groups(path, groups)
            store.commit(assignments, groups)
            self.assertEqual(store.initialize(assignments), {})
            self.assertEqual(training_feedback(path), {})
            groups[1] = 'good'
            self.save_groups(path, groups)
            store.commit(assignments, groups, reviewed_ids=[1])
            self.assertEqual(training_feedback(path), {1: 'good'})
            store.commit(assignments, groups, complete=True)
            self.assertEqual(len(training_feedback(path)), 3)
            # Structural edits invalidate only the affected cluster signatures.
            assignments[assignments == 1] = 0
            np.save(path / 'spike_clusters.npy', assignments)
            self.save_groups(path, {0: 'good', 2: 'noise'})
            store.commit(assignments, {0: 'good', 2: 'noise'})
            self.assertEqual(training_feedback(path), {2: 'noise'})
            # Edits outside the recorded save must not silently enter training.
            assignments[0] = 5
            np.save(path / 'spike_clusters.npy', assignments)
            self.assertEqual(training_feedback(path), {})

    def test_existing_curation_and_stale_predictions_are_protected(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            assignments = self.setup_dataset(path)
            self.save_groups(path, {0: 'noise', 1: 'good', 2: 'noise'})
            store = FeedbackStore(path)
            self.assertEqual(store.initialize(assignments), {})
            self.assertEqual(store.load()['initialization_policy'], 'existing_curation_preserved')
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            assignments = self.setup_dataset(path)
            assignments[0] = 1
            with self.assertRaises(ValueError):
                FeedbackStore(path).initialize(assignments)

    def test_unsaved_changes_cannot_be_confirmed(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            assignments = self.setup_dataset(path)
            store = FeedbackStore(path)
            store.initialize(assignments)
            store.initialized()
            with self.assertRaises(ValueError):
                store.commit(assignments, {0: 'good', 1: 'noise', 2: 'noise'}, complete=True)


if __name__ == '__main__':
    unittest.main()
