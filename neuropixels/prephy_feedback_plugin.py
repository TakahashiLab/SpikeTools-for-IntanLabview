"""Initialize prephy labels once and record saved human curation."""
import importlib.util
from pathlib import Path
import traceback

from phy import IPlugin, connect

spec = importlib.util.spec_from_file_location('_prephy_feedback_runtime', Path(__file__).with_name('prephy_feedback.py'))
feedback = importlib.util.module_from_spec(spec)
spec.loader.exec_module(feedback)


class PrephyFeedbackPlugin(IPlugin):
    def attach_to_controller(self, controller):
        @connect(sender=controller)
        def on_gui_ready(sender, gui):
            store = feedback.FeedbackStore(controller.model.dir_path)
            if not (store.directory / 'prephy_good_predictions.tsv').is_file() and store.load() is None:
                return
            supervisor = controller.supervisor
            touched, actions = set(), []
            pending_complete = [False]
            try:
                plan = store.initialize(supervisor.clustering.spike_clusters)
                for group in ('good', 'noise'):
                    ids = [cid for cid, label in plan.items() if label == group]
                    if ids:
                        supervisor.move(group, ids)
                store.initialized()
                if plan:
                    gui.status_message = 'Prephy: assigned predictions (good / otherwise noise); save with Ctrl+S.'
            except Exception:
                traceback.print_exc()
                gui.status_message = 'Prephy initialization failed; see console. Existing curation was preserved.'
                return

            @connect(sender=supervisor)
            def on_cluster(sender, up):
                description = str(getattr(up, 'description', ''))
                ids = list(getattr(up, 'metadata_changed', None) or [])
                history = getattr(up, 'history', None)
                if description == 'metadata_group':
                    if history == 'undo':
                        touched.difference_update(int(cid) for cid in ids)
                    else:
                        touched.update(int(cid) for cid in ids)
                pending_complete[0] = False
                actions.append({'description': description, 'history': history,
                                'cluster_ids': [int(cid) for cid in ids],
                                'value': getattr(up, 'metadata_value', None),
                                'added': [int(cid) for cid in (getattr(up, 'added', None) or [])],
                                'deleted': [int(cid) for cid in (getattr(up, 'deleted', None) or [])]})

            @connect(sender=supervisor, last=True)
            def on_save_clustering(sender, spike_clusters, groups, *labels):
                try:
                    store.commit(spike_clusters, groups, touched, pending_complete[0], actions)
                    touched.clear()
                    actions.clear()
                    pending_complete[0] = False
                    state = store.load()
                    gui.status_message = 'Prephy feedback saved: %d explicitly reviewed clusters.' % len(state['reviewed'])
                except Exception:
                    traceback.print_exc()
                    gui.status_message = 'Prephy feedback save failed; see console.'

            @gui.view_actions.add(alias='prephyreview', menu='Tools', submenu='Prephy feedback', shortcut='alt+shift+v')
            def confirm_selected_prephy_clusters():
                """Confirm selected clusters even when their predicted labels were correct."""
                touched.update(int(cid) for cid in supervisor.selected)
                supervisor.save()

            @gui.view_actions.add(alias='prephycomplete', menu='Tools', submenu='Prephy feedback', shortcut='alt+shift+c')
            def confirm_prephy_dataset_complete():
                """Confirm all current clusters have been manually checked, then save."""
                pending_complete[0] = True
                supervisor.save()
