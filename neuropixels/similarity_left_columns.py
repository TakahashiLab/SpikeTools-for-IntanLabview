"""Keep prephy metadata beside id in both phy cluster tables."""
from phy import IPlugin


def prephy_columns_first(columns):
    columns = list(dict.fromkeys(columns))
    by_name = {str(column).lower(): column for column in columns}
    priority = ('prephy_good_score', 'prephy_good_prediction',
                'prephy_feature_fraction', 'prephy_existing_group')
    prephy = [by_name[name] for name in priority if name in by_name]
    prephy.extend(column for column in columns
                  if str(column).lower().startswith('prephy') and column not in prephy)
    front = ['id'] + prephy
    if 'similarity' in columns:
        front.append('similarity')
    return front + [column for column in columns if column not in front]


class SimilarityLeftColumnsPlugin(IPlugin):
    def attach_to_controller(self, controller):
        from phy.cluster.supervisor import ClusterView

        if getattr(ClusterView, '_prephy_columns_patched', False):
            return
        original_reset = ClusterView._reset_table
        original_init = ClusterView.__init__

        def patched_fit_columns(self):
            from phy.gui.qt import QHeaderView
            header = self.table_view.horizontalHeader()
            # Continuous ResizeToContents and resizing every row are expensive
            # when SimilarityView replaces its data after each cluster choice.
            header.setSectionResizeMode(QHeaderView.Interactive)
            header.setResizeContentsPrecision(50)
            vertical = self.table_view.verticalHeader()
            vertical.setSectionResizeMode(QHeaderView.Fixed)
            vertical.setDefaultSectionSize(max(22, self.table_view.fontMetrics().height() + 8))
            schema = tuple(self.columns)
            if getattr(self, '_prephy_sized_columns', None) != schema:
                self.table_view.resizeColumnsToContents()
                self._prephy_sized_columns = schema

        def patched_reset(self, data=None, columns=(), sort=None):
            # Reorder only when phy itself creates/resets a table. Do not trigger
            # extra GUI resets, change sorting, or touch selection/merge actions.
            return original_reset(self, data=data,
                                  columns=prephy_columns_first(columns), sort=sort)

        def patched_init(self, *args, **kwargs):
            original_init(self, *args, **kwargs)
            # The local native Table implementation sets QTableView to NoFocus.
            # A row click then leaves keyboard focus in its QLineEdit filter,
            # which consumes G/K, Ctrl+Z and other standard phy shortcuts.
            from phy.gui.qt import Qt
            self.table_view.setFocusPolicy(Qt.StrongFocus)
            self.setFocusProxy(self.table_view)
            self.filter_edit.returnPressed.connect(self.table_view.setFocus)
            # Preserve coalescing while avoiding the previous 500 ms pause.
            self.debouncer.delay = 75

        def visible_sibling(self, row_id=None, direction='next'):
            # Predicted noise needs human review too. Follow the displayed
            # order, including masked/noise rows, without changing their group
            # or the styling mask. Explicit table filters still apply.
            selected = self.get_selected_ids()
            if row_id is None:
                row_id = selected[0] if selected else None
            visible = self._visible_ids()
            if row_id not in visible:
                return None
            index = visible.index(row_id) + (1 if direction == 'next' else -1)
            return visible[index] if 0 <= index < len(visible) else None

        def visible_first_or_last(self, which):
            visible = self._visible_ids()
            if visible:
                return self.select([visible[0] if which == 'first' else visible[-1]])
            return None

        # SimilarityView inherits this method, so both tables use the same order.
        ClusterView._reset_table = patched_reset
        ClusterView.__init__ = patched_init
        ClusterView._fit_columns = patched_fit_columns
        ClusterView.get_sibling_id = visible_sibling
        ClusterView._select_first_or_last = visible_first_or_last
        ClusterView._prephy_columns_patched = True
