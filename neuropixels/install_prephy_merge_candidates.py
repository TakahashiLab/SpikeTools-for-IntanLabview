"""Install merge suggestions while leaving the old learning plugin disabled."""
import ast
import argparse
from datetime import datetime
from pathlib import Path
import shutil


def install(phy_dir):
    source = Path(__file__).with_name('prephy_merge_candidates_plugin.py')
    destination = phy_dir / 'plugins' / source.name
    config = phy_dir / 'phy_config.py'
    original = config.read_text(encoding='utf-8')
    ast.parse(source.read_text(encoding='utf-8'))
    tree = ast.parse(original)
    assignment = next(node for node in tree.body if isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Attribute) and t.attr == 'plugins'
                and isinstance(t.value, ast.Attribute) and t.value.attr == 'TemplateGUI'
                for t in node.targets))
    names = ast.literal_eval(assignment.value)
    if 'PhyGoodPolicyPlugin' in names:
        raise RuntimeError('Disable PhyGoodPolicyPlugin first to avoid duplicate actions.')
    updated = original
    if 'PrephyMergeCandidatesPlugin' not in names:
        lines = original.splitlines(keepends=True)
        end = assignment.end_lineno - 1
        if lines[end].strip() != ']':
            raise RuntimeError('Expected multiline TemplateGUI plugin list.')
        lines.insert(end, "    'PrephyMergeCandidatesPlugin',\n")
        updated = ''.join(lines)
    ast.parse(updated)
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        shutil.copyfile(destination, destination.with_name(destination.name + '.' + stamp + '.bak'))
    shutil.copyfile(source, destination)
    if updated != original:
        shutil.copyfile(config, config.with_name(config.name + '.before_merge_candidates_' + stamp + '.bak'))
        config.write_text(updated, encoding='utf-8')
    print('Enabled PrephyMergeCandidatesPlugin. Restart Phy to use merge suggestions.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--phy-dir', type=Path, default=Path.home() / '.phy')
    install(parser.parse_args().phy_dir)
