"""Test real Qt keyboard delivery offscreen, before and after the plugin fix."""
import importlib.util
import os
import sys
from pathlib import Path
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

from PyQt5.QtCore import Qt
from PyQt5.QtTest import QTest
from PyQt5.QtWidgets import QApplication, QAction, QMainWindow
from phy.cluster.supervisor import ClusterView, SimilarityView

app = QApplication.instance() or QApplication([])
app.setQuitOnLastWindowClosed(False)


def exercise(view_type, fixed):
    window = QMainWindow()
    data = [{'id': cid, 'n_spikes': 10, 'prephy_good_score': cid / 10, 'similarity': 0.8,
             'group': 'noise', 'is_masked': True} for cid in (1, 2, 3)]
    view = view_type(window, data=data, columns=['id', 'n_spikes', 'prephy_good_score', 'similarity'], sort=('id', 'asc'))
    window.setCentralWidget(view)
    events = []
    def add(key, name, callback=None):
        action = QAction(name, window)
        action.setShortcut(key)
        def run():
            events.append(name)
            if callback:
                callback()
        action.triggered.connect(run)
        window.addAction(action)
    add('Down', 'down', view.next)
    add('Up', 'up', view.previous)
    add('g', 'merge')
    add('k', 'split')
    add('Ctrl+S', 'save')
    add('Ctrl+Z', 'undo')
    add('Alt+G', 'good')
    window.resize(700, 300)
    window.show()
    window.activateWindow()
    QTest.qWait(100)
    view.filter_edit.setFocus()
    app.processEvents()
    index = view._proxy.index(0, 0)
    QTest.mouseClick(view.table_view.viewport(), Qt.LeftButton, pos=view.table_view.visualRect(index).center())
    app.processEvents()
    assert view.get_selected_ids() == [1]
    focus_after_click = app.focusWidget()
    QTest.keyClick(focus_after_click, Qt.Key_G)
    app.processEvents()
    if not fixed:
        assert focus_after_click is view.filter_edit, 'Expected to reproduce the NoFocus filter trap.'
        assert 'merge' not in events and view.filter_edit.text() == 'g'
        print('REPRODUCED: row click leaves focus in filter; G is typed instead of merge.', flush=True)
    else:
        assert focus_after_click is view.table_view
        QTest.keyClick(app.focusWidget(), Qt.Key_Down)
        QTest.qWait(150)
        assert view.get_selected_ids() == [2], (events, view.get_selected_ids())
        QTest.keyClick(app.focusWidget(), Qt.Key_Up)
        QTest.qWait(150)
        assert view.get_selected_ids() == [1]
        QTest.keyClick(app.focusWidget(), Qt.Key_K)
        QTest.keyClick(app.focusWidget(), Qt.Key_S, Qt.ControlModifier)
        QTest.keyClick(app.focusWidget(), Qt.Key_Z, Qt.ControlModifier)
        QTest.keyClick(app.focusWidget(), Qt.Key_G, Qt.AltModifier)
        assert events == ['merge', 'down', 'up', 'split', 'save', 'undo', 'good'], events
        assert view.filter_edit.text() == ''
        assert view.columns[:2] == ['id', 'prephy_good_score']
        # All-noise data must support first/last and both arrow directions.
        view.select([])
        view.first()
        assert view.get_selected_ids() == [1]
        view.last()
        assert view.get_selected_ids() == [3]
        view.next()
        assert view.get_selected_ids() == [3]  # No wrap past the boundary.
        view.sort_by('id', 'desc')
        view.next()
        assert view.get_selected_ids() == [2]
        view.filter('id == 2')
        view.first()
        assert view.get_selected_ids() == [2]
        assert view.get_sibling_id(2) is None
        view.filter('')
        # Navigation never changes the noise groups or visual masks.
        assert all(row['group'] == 'noise' and row['is_masked'] for row in view._model._rows)
        # Intentional filtering still works; Enter returns keyboard focus to rows.
        view.filter_edit.setFocus()
        QTest.keyClick(view.filter_edit, Qt.Key_Return)
        app.processEvents()
        assert app.focusWidget() is view.table_view
        print('PASS: ' + view_type.__name__ + ' noise arrows, first/last, sort/filter/boundaries, G/K, Ctrl+S/Z, Alt+G and prephy columns.', flush=True)
    window.close()
    app.processEvents()


mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
if mode in ('all', 'baseline'):
    exercise(ClusterView, fixed=False)
spec = importlib.util.spec_from_file_location('prephy_columns_test', Path(__file__).with_name('similarity_left_columns.py'))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.SimilarityLeftColumnsPlugin().attach_to_controller(None)
if mode in ('all', 'cluster'):
    exercise(ClusterView, fixed=True)
if mode in ('all', 'similarity'):
    exercise(SimilarityView, fixed=True)
