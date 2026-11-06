"""Manual good/noise buffers only; no training, predictions or auto-assignment."""
import traceback
from phy import IPlugin


class LabelBuffers:
    def __init__(self):
        self.good = set()
        self.noise = set()

    def prune(self, active):
        active = set(int(cid) for cid in active)
        self.good.intersection_update(active)
        self.noise.intersection_update(active)

    def add(self, group, selected, active):
        self.prune(active)
        active = set(int(cid) for cid in active)
        selected = active.intersection(int(cid) for cid in selected)
        target = getattr(self, group)
        other = self.noise if group == 'good' else self.good
        before = len(target)
        target.update(selected)
        other.difference_update(selected)
        return len(target) - before

    def remove(self, group, selected):
        target = getattr(self, group)
        before = len(target)
        target.difference_update(int(cid) for cid in selected)
        return before - len(target)

    def apply(self, group, supervisor):
        self.prune(supervisor.clustering.cluster_ids)
        target = getattr(self, group)
        candidates = sorted(cid for cid in target if supervisor.cluster_meta.get('group', cid) != group)
        if candidates:
            # One native phy action: supports undo and emits ordinary manual
            # group events consumed by PrephyFeedbackPlugin when saved.
            supervisor.move(group, candidates)
        target.clear()
        return len(candidates)


class PrephyLabelBufferPlugin(IPlugin):
    def attach_to_controller(self, controller):
        from phy import connect

        @connect(sender=controller)
        def on_gui_ready(sender, gui):
            supervisor = controller.supervisor
            buffers = LabelBuffers()
            self.buffers = buffers

            def perform(operation, group):
                try:
                    selected = list(supervisor.selected)
                    if operation == 'add':
                        count = buffers.add(group, selected, supervisor.clustering.cluster_ids)
                    elif operation == 'remove':
                        buffers.prune(supervisor.clustering.cluster_ids)
                        count = buffers.remove(group, selected)
                    else:
                        count = buffers.apply(group, supervisor)
                    gui.status_message = ('Label buffer: %s %s=%d; good=%d noise=%d%s' %
                                          (operation, group, count, len(buffers.good), len(buffers.noise),
                                           '; Ctrl+S to save.' if operation == 'apply' else ''))
                except Exception:
                    traceback.print_exc()
                    gui.status_message = 'Label buffer operation failed; see console.'

            @gui.view_actions.add(alias='addgood', menu='Tools', submenu='Label buffers', shortcut='alt+shift+j', toolbar=True)
            def buffer_selected_good():
                """Add selected clusters to the good buffer without changing labels."""
                perform('add', 'good')

            @gui.view_actions.add(alias='unbuffergood', menu='Tools', submenu='Label buffers', shortcut='alt+shift+k', toolbar=True)
            def unbuffer_selected_good():
                """Remove selected clusters from the good buffer."""
                perform('remove', 'good')

            @gui.view_actions.add(alias='applygood', menu='Tools', submenu='Label buffers', shortcut='alt+shift+l', toolbar=True)
            def apply_good_buffer():
                """Assign all buffered good labels in one undoable action."""
                perform('apply', 'good')

            @gui.view_actions.add(alias='addnoise', menu='Tools', submenu='Label buffers', shortcut='alt+shift+a', toolbar=True)
            def buffer_selected_noise():
                """Add selected clusters to the noise buffer without changing labels."""
                perform('add', 'noise')

            @gui.view_actions.add(alias='unbuffernoise', menu='Tools', submenu='Label buffers', shortcut='alt+shift+u', toolbar=True)
            def unbuffer_selected_noise():
                """Remove selected clusters from the noise buffer."""
                perform('remove', 'noise')

            @gui.view_actions.add(alias='applynoise', menu='Tools', submenu='Label buffers', shortcut='alt+shift+b', toolbar=True)
            def apply_noise_buffer():
                """Assign all buffered noise labels in one undoable action."""
                perform('apply', 'noise')
