"""Install only the feedback plugin and its dependency; retain disabled policy."""
import ast
from datetime import datetime
from pathlib import Path
import shutil

source = Path(__file__).resolve().parent
destination = Path.home() / '.phy' / 'plugins'
config = destination.parent / 'phy_config.py'
destination.mkdir(parents=True, exist_ok=True)
stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
for name in ('prephy_feedback.py', 'prephy_feedback_plugin.py'):
    ast.parse((source / name).read_text(encoding='utf-8'))
    target = destination / name
    if target.exists():
        shutil.copyfile(target, target.with_name(target.name + '.' + stamp + '.bak'))
    shutil.copyfile(source / name, target)
original = config.read_text(encoding='utf-8')
tree = ast.parse(original)
assignment = next(node for node in tree.body if isinstance(node, ast.Assign)
                  and any(isinstance(target, ast.Attribute) and target.attr == 'plugins'
                          and isinstance(target.value, ast.Attribute) and target.value.attr == 'TemplateGUI'
                          for target in node.targets))
names = ast.literal_eval(assignment.value)
if 'PrephyFeedbackPlugin' not in names:
    lines = original.splitlines(keepends=True)
    end = assignment.end_lineno - 1
    # Insert before the closing bracket without re-enabling disabled plugins.
    if lines[end].strip() != ']':
        raise RuntimeError('Expected a multiline TemplateGUI plugin list.')
    lines.insert(end, "    'PrephyFeedbackPlugin',\n")
    updated = ''.join(lines)
    ast.parse(updated)
    shutil.copyfile(config, config.with_name(config.name + '.before_feedback_' + stamp + '.bak'))
    config.write_text(updated, encoding='utf-8')
print('Installed PrephyFeedbackPlugin; PhyGoodPolicyPlugin remains disabled.')
