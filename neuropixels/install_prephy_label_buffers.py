import ast
from datetime import datetime
from pathlib import Path
import shutil

source = Path(__file__).with_name('prephy_label_buffer_plugin.py')
target = Path.home() / '.phy' / 'plugins' / 'prephy_label_buffer_plugin.py'
config = target.parent.parent / 'phy_config.py'
target.parent.mkdir(parents=True, exist_ok=True)
stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
ast.parse(source.read_text(encoding='utf-8'))
original = config.read_text(encoding='utf-8')
tree = ast.parse(original)
assignment = next(node for node in tree.body if isinstance(node, ast.Assign)
                  and any(isinstance(item, ast.Attribute) and item.attr == 'plugins'
                          and isinstance(item.value, ast.Attribute) and item.value.attr == 'TemplateGUI'
                          for item in node.targets))
names = ast.literal_eval(assignment.value)
assert 'PhyGoodPolicyPlugin' not in names, 'Old learning plugin must stay disabled.'
updated = original
if 'PrephyLabelBufferPlugin' not in names:
    lines = original.splitlines(keepends=True)
    end = assignment.end_lineno - 1
    assert lines[end].strip() == ']'
    lines.insert(end, "    'PrephyLabelBufferPlugin',\n")
    updated = ''.join(lines)
ast.parse(updated)
if target.exists():
    shutil.copyfile(target, target.with_name(target.name + '.' + stamp + '.bak'))
shutil.copyfile(source, target)
if updated != original:
    shutil.copyfile(config, config.with_name(config.name + '.before_label_buffers_' + stamp + '.bak'))
    config.write_text(updated, encoding='utf-8')
print('Enabled six manual buffer actions only. Old learning plugin remains disabled.')
