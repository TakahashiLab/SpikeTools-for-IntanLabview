"""Manual merge suggestions only; no model loading, training or automatic labels."""
import json
import traceback
from pathlib import Path
import numpy as np
from phy import IPlugin


class MergeCandidates:
    def __init__(self, controller, similarity=0.81, refractory_ms=1.2):
        self.controller = controller
        self.similarity = float(similarity)
        self.refractory_s = float(refractory_ms) / 1000
        if not np.isfinite(self.similarity) or not 0 <= self.similarity <= 1:
            raise ValueError('Invalid merge similarity threshold')
        if not np.isfinite(self.refractory_s) or self.refractory_s < 0:
            raise ValueError('Invalid refractory interval')
        self.rejected = set()

    def selected_pair(self):
        supervisor = self.controller.supervisor
        selected = list(dict.fromkeys(int(c) for c in supervisor.selected))
        active = set(int(c) for c in supervisor.clustering.cluster_ids)
        if len(selected) == 2 and set(selected) <= active:
            return tuple(sorted(selected))
        return None

    def ranked(self, reference=None):
        controller = self.controller
        supervisor = controller.supervisor
        active = sorted(int(c) for c in supervisor.clustering.cluster_ids
                        if supervisor.cluster_meta.get('group', int(c)) != 'noise')
        if reference is not None and reference not in active:
            return []
        matrix = np.asarray(controller.model.similar_templates)
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            raise ValueError('similar_templates must be a square matrix')
        # Cluster IDs change after merge/split. Use the current template membership,
        # never index the original similarity matrix by a curated cluster ID.
        templates = {c: np.flatnonzero(controller.get_template_counts(c)) for c in active}
        times = {}

        def spike_times(c):
            if c not in times:
                value = controller.get_spike_times(c)
                times[c] = np.sort(np.asarray(getattr(value, 'data', value), dtype=float))
            return times[c]

        pairs = []
        lefts = active if reference is None else [reference]
        for left in lefts:
            for right in active:
                if right == left or (reference is None and right < left):
                    continue
                pair = tuple(sorted((left, right)))
                if pair in self.rejected:
                    continue
                a, b = templates[left], templates[right]
                if not len(a) or not len(b):
                    continue
                score = float(np.max(matrix[np.ix_(a, b)]))
                if not np.isfinite(score) or score < self.similarity:
                    continue
                pairs.append((pair, score))
        pairs.sort(key=lambda item: (-item[1], item[0]))
        result = []
        for pair, score in pairs:
            a, b = spike_times(pair[0]), spike_times(pair[1])
            lo = np.searchsorted(b, a - self.refractory_s, side='left')
            hi = np.searchsorted(b, a + self.refractory_s, side='right')
            if not np.any(hi > lo):
                result.append((pair, score))
        return result


class PrephyMergeCandidatesPlugin(IPlugin):
    def attach_to_controller(self, controller):
        from phy import connect
        config_path = Path.home() / '.phy' / 'good_policy_config.json'
        config = json.loads(config_path.read_text(encoding='utf-8-sig')) if config_path.exists() else {}
        self.candidates = MergeCandidates(controller,
            config.get('merge_similarity_threshold', 0.81), config.get('xcorr_refractory_ms', 1.2))

        @connect(sender=controller)
        def on_gui_ready(sender, gui):
            engine = self.candidates
            supervisor = controller.supervisor

            def perform(action):
                try:
                    action()
                except Exception:
                    traceback.print_exc()
                    gui.status_message = 'Merge candidates failed; see console.'

            def next_pair(prefix=''):
                ranked = engine.ranked()
                if not ranked:
                    gui.status_message = prefix + 'No remaining merge candidates (similarity >= %.3f, refractory %.2f ms).' % (engine.similarity, engine.refractory_s * 1000)
                    return
                pair, score = ranked[0]
                supervisor.select(list(pair))
                gui.status_message = prefix + 'Review %s similarity=%.4f; :domerge or :rejectmerge' % (pair, score)

            @gui.view_actions.add(alias='nextmerge', menu='Tools', submenu='Merge candidates', shortcut='alt+shift+m', toolbar=True)
            def select_next_merge_pair():
                """Select the highest similarity eligible pair for manual review."""
                perform(next_pair)

            @gui.view_actions.add(alias='mergesuggest', menu='Tools', submenu='Merge candidates')
            def select_merge_candidates(n='20'):
                """Select up to N candidates for the first selected cluster."""
                def run():
                    selected = list(supervisor.selected)
                    if not selected:
                        gui.status_message = 'Select a reference cluster before :mergesuggest'
                        return
                    reference = int(selected[0])
                    ranked = engine.ranked(reference)[:max(0, int(n))]
                    others = [pair[1] if pair[0] == reference else pair[0] for pair, _ in ranked]
                    supervisor.select([reference] + others)
                    gui.status_message = 'Selected %d candidates for %d.' % (len(others), reference)
                perform(run)

            @gui.view_actions.add(alias='domerge', menu='Tools', submenu='Merge candidates', shortcut='alt+shift+return', toolbar=True)
            def merge_selected_pair():
                """Merge exactly the two currently selected clusters, then suggest next."""
                def run():
                    pair = engine.selected_pair()
                    if pair is None:
                        gui.status_message = 'Select exactly two current clusters before :domerge'
                        return
                    supervisor.merge(list(pair))
                    next_pair('Merged %s; Ctrl+S to save. ' % (pair,))
                perform(run)

            @gui.view_actions.add(alias='rejectmerge', menu='Tools', submenu='Merge candidates', shortcut='ctrl+shift+r', toolbar=True)
            def reject_selected_pair():
                """Skip the selected pair for this session and suggest the next pair."""
                def run():
                    pair = engine.selected_pair()
                    if pair is None:
                        gui.status_message = 'Select exactly two current clusters before :rejectmerge'
                        return
                    engine.rejected.add(pair)
                    next_pair('Rejected %s. ' % (pair,))
                perform(run)
