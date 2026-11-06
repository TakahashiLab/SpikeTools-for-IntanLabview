import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest
from phy import emit, unconnect

spec = importlib.util.spec_from_file_location('prephy_label_buffers_test', Path(__file__).with_name('prephy_label_buffer_plugin.py'))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class BufferTests(unittest.TestCase):
    def test_manual_buffer_workflow_without_predictions_or_training(self):
        groups = {1: 'noise', 2: 'noise', 3: 'good'}
        calls = []
        supervisor = SimpleNamespace(selected=[1, 2], clustering=SimpleNamespace(cluster_ids=[1, 2, 3]),
                                     cluster_meta=SimpleNamespace(get=lambda field, cid: groups[cid]))
        def move(group, ids):
            calls.append((group, list(ids)))
            for cid in ids:
                groups[cid] = group
        supervisor.move = move
        callbacks, keys = {}, {}
        class Actions:
            def add(self, **kwargs):
                def decorate(callback):
                    callbacks[kwargs['alias']] = callback
                    keys[kwargs['alias']] = kwargs['shortcut']
                    return callback
                return decorate
        controller = SimpleNamespace(supervisor=supervisor)
        gui = SimpleNamespace(view_actions=Actions(), status_message='')
        plugin = module.PrephyLabelBufferPlugin()
        plugin.attach_to_controller(controller)
        try:
            emit('gui_ready', controller, gui)
            self.assertEqual(len(callbacks), 6)
            self.assertEqual(len(set(keys.values())), 6)
            callbacks['addgood']()
            self.assertEqual(groups, {1: 'noise', 2: 'noise', 3: 'good'})
            self.assertEqual(plugin.buffers.good, {1, 2})
            supervisor.selected = [2]
            callbacks['unbuffergood']()
            self.assertEqual(plugin.buffers.good, {1})
            callbacks['addgood']()
            callbacks['addnoise']()  # Opposite buffer wins, without assignment.
            self.assertEqual(plugin.buffers.good, {1})
            self.assertEqual(plugin.buffers.noise, {2})
            callbacks['unbuffernoise']()
            self.assertEqual(plugin.buffers.noise, set())
            callbacks['applygood']()
            self.assertEqual(calls, [('good', [1])])
            self.assertEqual(plugin.buffers.good, set())
            supervisor.selected = [1, 3]
            callbacks['addnoise']()
            supervisor.clustering.cluster_ids = [1, 2]  # Cluster 3 merged/deleted.
            callbacks['applynoise']()
            self.assertEqual(calls[-1], ('noise', [1]))
            self.assertEqual(plugin.buffers.noise, set())
            callbacks['applynoise']()
            self.assertEqual(len(calls), 2)
        finally:
            unconnect(controller)

    def test_failure_retains_pending_labels(self):
        buffers = module.LabelBuffers()
        buffers.add('good', [1], [1])
        def fail(*args):
            raise RuntimeError('move failed')
        supervisor = SimpleNamespace(clustering=SimpleNamespace(cluster_ids=[1]),
                                     cluster_meta=SimpleNamespace(get=lambda *args: 'noise'), move=fail)
        with self.assertRaises(RuntimeError):
            buffers.apply('good', supervisor)
        self.assertEqual(buffers.good, {1})


if __name__ == '__main__':
    unittest.main()
