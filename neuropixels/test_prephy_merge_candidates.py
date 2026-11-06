import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace as NS
import numpy as np
from prephy_merge_candidates_plugin import MergeCandidates
from install_prephy_merge_candidates import install


class MergeCandidatesTest(unittest.TestCase):
    def setUp(self):
        self.groups = {10: 'good', 20: 'mua', 30: 'noise', 40: 'good'}
        self.counts = {10: [1, 0, 0], 20: [0, 1, 0], 30: [0, 0, 1], 40: [1, 1, 0]}
        self.times = {10: [0., 1.], 20: [.1, 1.1], 30: [.2], 40: [0., .1, 1., 1.1]}
        self.sup = NS(selected=[10, 20], clustering=NS(cluster_ids=[10, 20, 30]),
                      cluster_meta=NS(get=lambda key, c: self.groups[c]))
        self.controller = NS(supervisor=self.sup,
            model=NS(similar_templates=np.array([[1., .9, .99], [.9, 1., .95], [.99, .95, 1.]])),
            get_template_counts=lambda c: self.counts[c],
            get_spike_times=lambda c: NS(data=np.array(self.times[c])))
        self.engine = MergeCandidates(self.controller)

    def test_current_ids_and_noise(self):
        self.assertEqual(self.engine.ranked(), [((10, 20), .9)])
        self.assertEqual(self.engine.ranked(20), [((10, 20), .9)])
        self.assertEqual(self.engine.ranked(30), [])

    def test_refractory_and_reject(self):
        self.times[20] = [.0012]
        self.assertEqual(self.engine.ranked(), [])
        self.times[20] = [.1]
        self.engine.rejected.add((10, 20))
        self.assertEqual(self.engine.ranked(), [])

    def test_membership_after_merge_or_split(self):
        self.sup.clustering.cluster_ids = [40, 20]
        self.times[40] = [.3]
        self.assertEqual(self.engine.ranked(), [((20, 40), 1.)])
        self.sup.selected = [10, 20]
        self.assertIsNone(self.engine.selected_pair())
        self.sup.selected = [40, 20]
        self.assertEqual(self.engine.selected_pair(), (20, 40))
        self.sup.selected = [40]
        self.assertIsNone(self.engine.selected_pair())

    def test_threshold_and_nan(self):
        self.engine.similarity = .95
        self.assertEqual(self.engine.ranked(), [])
        self.controller.model.similar_templates[0, 1] = np.nan
        self.assertEqual(self.engine.ranked(), [])

    def test_gui_actions_require_explicit_merge(self):
        from phy import emit
        from unittest.mock import patch
        from prephy_merge_candidates_plugin import PrephyMergeCandidatesPlugin
        actions = {}
        def add(**options):
            def register(fn):
                actions[options['alias']] = fn
                return fn
            return register
        gui = NS(view_actions=NS(add=add), status_message='')
        merged = []
        self.sup.select = lambda ids: setattr(self.sup, 'selected', ids)
        def merge(ids):
            merged.append(ids)
            self.sup.clustering.cluster_ids = [40, 30]
            self.sup.selected = [40]
        self.sup.merge = merge
        plugin = PrephyMergeCandidatesPlugin()
        with tempfile.TemporaryDirectory() as tmp, patch('prephy_merge_candidates_plugin.Path.home', return_value=Path(tmp)):
            plugin.attach_to_controller(self.controller)
        emit('gui_ready', self.controller, gui)
        self.assertEqual(set(actions), {'nextmerge', 'mergesuggest', 'domerge', 'rejectmerge'})
        actions['nextmerge']()
        self.assertEqual(self.sup.selected, [10, 20])
        self.assertEqual(merged, [])
        actions['rejectmerge']()
        self.assertEqual(merged, [])
        self.assertEqual(plugin.candidates.ranked(), [])
        plugin.candidates.rejected.clear()
        actions['mergesuggest']('1')
        self.assertEqual(self.sup.selected, [10, 20])
        self.sup.selected = [10]
        actions['domerge']()
        self.assertEqual(merged, [])
        self.sup.selected = [10, 20]
        actions['domerge']()
        self.assertEqual(merged, [[10, 20]])

    def test_installer_preserves_plugins_and_backup(self):
        with tempfile.TemporaryDirectory() as tmp:
            phy_dir = Path(tmp)
            config = phy_dir / 'phy_config.py'
            text = "c = get_config()\nc.TemplateGUI.plugins = [\n    'PrephyFeedbackPlugin',\n    'PrephyLabelBufferPlugin',\n]\n"
            config.write_text(text, encoding='utf-8')
            install(phy_dir)
            saved = config.read_text(encoding='utf-8')
            self.assertIn("'PrephyMergeCandidatesPlugin'", saved)
            self.assertIn("'PrephyFeedbackPlugin'", saved)
            self.assertNotIn("'PhyGoodPolicyPlugin'", saved)
            self.assertEqual(next(phy_dir.glob('*.bak')).read_text(encoding='utf-8'), text)
            install(phy_dir)
            self.assertEqual(config.read_text(encoding='utf-8'), saved)


if __name__ == '__main__':
    unittest.main()
