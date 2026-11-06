"""Exercise real phy events without opening a GUI or touching real data."""
import importlib.util
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
from phy import emit, unconnect
from prephy_feedback import training_feedback, read_groups
from test_prephy_feedback import FeedbackTests
from test_prephy_label_buffers import module as buffer_module

spec = importlib.util.spec_from_file_location('prephy_feedback_plugin', Path(__file__).with_name('prephy_feedback_plugin.py'))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class PluginEventsTest(unittest.TestCase):
    def test_initialization_save_and_explicit_review(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            assignments = FeedbackTests().setup_dataset(path)
            class Supervisor:
                clustering = SimpleNamespace(spike_clusters=assignments)
                selected = [0]
                groups = {0: 'unsorted', 1: 'unsorted', 2: 'unsorted'}

                def move(self, group, ids):
                    for cid in ids:
                        self.groups[cid] = group
                    emit('cluster', self, SimpleNamespace(description='metadata_group', metadata_changed=ids,
                                                         metadata_value=group, history=None, added=[], deleted=[]))

                def save(self):
                    FeedbackTests().save_groups(path, self.groups)
                    np.save(path / 'spike_clusters.npy', self.clustering.spike_clusters)
                    emit('save_clustering', self, self.clustering.spike_clusters, self.groups)

            supervisor = Supervisor()
            supervisor.clustering.cluster_ids = np.unique(assignments)
            supervisor.cluster_meta = SimpleNamespace(get=lambda field, cid: supervisor.groups[cid])
            callbacks = {}
            class Actions:
                def add(self, **kwargs):
                    def decorate(callback):
                        callbacks[kwargs['alias']] = callback
                        return callback
                    return decorate
            controller = SimpleNamespace(model=SimpleNamespace(dir_path=str(path)), supervisor=supervisor)
            gui = SimpleNamespace(view_actions=Actions(), status_message='')
            module.PrephyFeedbackPlugin().attach_to_controller(controller)
            buffer_module.PrephyLabelBufferPlugin().attach_to_controller(controller)
            try:
                emit('gui_ready', controller, gui)
                self.assertEqual(supervisor.groups, {0: 'good', 1: 'noise', 2: 'noise'})
                supervisor.save()
                self.assertEqual(training_feedback(path), {})
                supervisor.selected = [1]
                callbacks['addgood']()
                self.assertEqual(supervisor.groups[1], 'noise')
                self.assertEqual(training_feedback(path), {})
                callbacks['applygood']()
                supervisor.save()
                self.assertEqual(training_feedback(path), {1: 'good'})
                supervisor.selected = [0]
                callbacks['prephyreview']()
                self.assertEqual(training_feedback(path), {0: 'good', 1: 'good'})
                callbacks['prephycomplete']()
                self.assertEqual(training_feedback(path), read_groups(path))
            finally:
                unconnect(controller, supervisor)


if __name__ == '__main__':
    unittest.main()
