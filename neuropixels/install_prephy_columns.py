"""Validate and install the existing column plugin with a backup."""
import ast
from datetime import datetime
from pathlib import Path
import sys
import types

source = Path(__file__).with_name('similarity_left_columns.py')
code = source.read_text(encoding='utf-8')
ast.parse(code)

# Exercise the installed phy method interface without opening a GUI or dataset.
phy = types.ModuleType('phy')
phy.IPlugin = object
supervisor = types.ModuleType('phy.cluster.supervisor')

class ClusterView:
    def _reset_table(self, data=None, columns=(), sort=None):
        return data, columns, sort

class SimilarityView(ClusterView):
    pass

supervisor.ClusterView = ClusterView
sys.modules['phy'] = phy
sys.modules['phy.cluster.supervisor'] = supervisor
namespace = {}
exec(compile(code, str(source), 'exec'), namespace)
plugin = namespace['SimilarityLeftColumnsPlugin']()
plugin.attach_to_controller(None)
original_patch = ClusterView._reset_table
plugin.attach_to_controller(None)
assert ClusterView._reset_table is original_patch
data = [{'id': 1}]
sort = ('prephy_good_score', 'desc')
for view in (object.__new__(ClusterView), object.__new__(SimilarityView)):
    returned, columns, sorting = view._reset_table(data, ['id', 'ContamPct', 'similarity', 'prephy_good_prediction', 'prephy_good_score', 'prephy_extra', 'id'], sort)
    assert columns == ['id', 'prephy_good_score', 'prephy_good_prediction', 'prephy_extra', 'similarity', 'ContamPct']
    assert returned is data and sorting is sort
    assert view._reset_table(columns=['id', 'amp'])[1] == ['id', 'amp']
print('PASS: both tables, duplicate columns, unchanged sorting/data, repeated attachment.')

if '--apply' in sys.argv:
    target = Path.home() / '.phy' / 'plugins' / 'similarity_left_columns.py'
    backup = target.with_name(target.name + '.before_prephy_columns_' + datetime.now().strftime('%Y%m%d_%H%M%S') + '.bak')
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        backup.write_bytes(target.read_bytes())
    target.write_text(code, encoding='utf-8')
    assert target.read_text(encoding='utf-8') == code
    print('Updated: ' + str(target))
    print('Backup: ' + str(backup))
