"""Train and apply a good-cluster classifier without starting phy.

Training labels come exclusively from cluster_group.tsv, never KSLabel.
Prediction writes separate metadata TSV files and never changes curation.
"""
import argparse
import ast
from collections import Counter
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import sys

import joblib
import numpy as np
from prephy_feedback import array_hash, training_feedback
from sklearn import __version__ as sklearn_version
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (average_precision_score, confusion_matrix,
                             precision_score, recall_score, roc_auc_score)
from sklearn.model_selection import GroupShuffleSplit

VERSION = '1.1.0'
FEATURES = (
    'n_spikes', 'firing_rate', 'amplitude', 'contam_pct', 'isi_viol', 'mergemax',
    'duration', 'halfwidth', 'pt_ratio', 'repolarization_slope', 'recovery_slope',
    'snr', 'spread', 'velocity_above', 'velocity_below',
)
ALIASES = {
    'n_spikes': ('n_spikes', 'n_spks', 'spike_count'),
    'firing_rate': ('firing_rate', 'fr'),
    'amplitude': ('Amplitude', 'amplitude', 'amp'),
    'contam_pct': ('ContamPct', 'contam_pct', 'contamination'),
    'isi_viol': ('isi_viol', 'isi_violations'),
    'mergemax': ('MergeMax', 'mergemax', 'merge_max'),
    'pt_ratio': ('PT_ratio', 'pt_ratio', 'peak_trough_ratio'),
    'halfwidth': ('halfwidth', 'half_width'),
}
SOURCE_FILES = ('cluster_info.tsv', 'cluster_metrics.tsv', 'metrics.csv',
                'waveform_metrics.csv', 'cluster_ContamPct.tsv',
                'cluster_Amplitude.tsv', 'cluster_isi_viol.tsv',
                'cluster_MergeMax.tsv', 'cluster_mergemax.tsv',
                'cluster_group.tsv', 'spike_clusters.npy', 'spike_times.npy',
                'similar_templates.npy', 'params.py')


def log(message):
    print(message, flush=True)


def read_rows(path):
    if not path.is_file():
        return []
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter=',' if path.suffix == '.csv' else '\t'))


def cluster_id(row):
    for key in ('cluster_id', 'id', 'cluster'):
        try:
            return int(row[key])
        except (KeyError, ValueError, TypeError):
            continue
    return None


def number(row, keys):
    for key in keys:
        try:
            value = float(row.get(key, ''))
            if np.isfinite(value):
                return value
        except (ValueError, TypeError):
            pass
    return np.nan


def manual_labels(directory):
    # Do not fall back to cluster_KSLabel.tsv or group columns in feature files.
    return {cluster_id(row): str(row.get('group', '')).strip().lower()
            for row in read_rows(directory / 'cluster_group.tsv')
            if cluster_id(row) is not None}


def discover(roots):
    """Only traverse run/probe directories; avoid huge .phy/joblib trees."""
    found = {}
    for root in roots:
        root = Path(root).expanduser().resolve()
        if not root.is_dir():
            raise ValueError('Input directory does not exist: ' + str(root))
        for current, dirs, files in os.walk(root):
            dirs[:] = sorted(d for d in dirs if d not in ('.phy', '__pycache__', 'recovery', 'tools', 'prephy_feedback')
                             and not d.endswith('_orig') and '_failed' not in d)
            directory = Path(current)
            if 'params.py' in files and 'spike_clusters.npy' in files:
                found[str(directory.resolve()).casefold()] = directory.resolve()
                dirs[:] = []
    if not found:
        raise ValueError('No phy datasets containing params.py and spike_clusters.npy found.')
    return sorted(found.values())


def recording_group(directory):
    for parent in directory.parents:
        if re.fullmatch(r'(?:catgt_)?[^\\/]+_g\d+', parent.name):
            return str(parent).casefold()
    return str(directory.parent).casefold()


def parameters(directory):
    values = {}
    # Parse literal assignments without executing an arbitrary params.py.
    tree = ast.parse((directory / 'params.py').read_text(encoding='utf-8-sig'))
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                values[node.targets[0].id] = ast.literal_eval(node.value)
            except (ValueError, TypeError):
                continue
    return values


def recording_duration(directory, params, spike_times):
    dat_path = params.get('dat_path')
    if isinstance(dat_path, (list, tuple)) and len(dat_path) == 1:
        dat_path = dat_path[0]
    if isinstance(dat_path, str):
        binary = Path(dat_path)
        if not binary.is_absolute():
            binary = directory / binary
        meta = binary.with_suffix('.meta')
        if meta.is_file():
            for line in meta.read_text(encoding='utf-8').splitlines():
                if line.startswith('fileTimeSecs='):
                    value = float(line.split('=', 1)[1])
                    if value > 0:
                        return value, 'metadata'
        if binary.is_file() and params.get('n_channels_dat'):
            sample_bytes = np.dtype(params.get('dtype', 'int16')).itemsize * int(params['n_channels_dat'])
            value = (binary.stat().st_size - int(params.get('offset', 0))) / sample_bytes / float(params['sample_rate'])
            if value > 0:
                return value, 'binary_size'
    if spike_times.size:
        return (float(np.max(spike_times)) + 1) / float(params['sample_rate']), 'last_spike'
    raise ValueError('Cannot determine duration: ' + str(directory))


def extract(directory):
    assignments = np.load(directory / 'spike_clusters.npy', mmap_mode='r').reshape(-1)
    if not assignments.size:
        raise ValueError('No active spikes: ' + str(directory))
    ids, counts = np.unique(assignments, return_counts=True)
    ids = ids.astype(np.int64)
    rows = {int(cid): {} for cid in ids}
    # Use the same sources for training and new, not-yet-opened-in-phy data.
    # Current spike membership takes precedence over stale cluster_info counts.
    for filename in ('cluster_info.tsv', 'cluster_metrics.tsv', 'metrics.csv', 'waveform_metrics.csv'):
        for row in read_rows(directory / filename):
            cid = cluster_id(row)
            if cid not in rows:
                continue
            for feature in FEATURES:
                value = number(row, ALIASES.get(feature, (feature,)))
                if np.isfinite(value):
                    rows[cid][feature] = value
    for feature, filenames in (
        ('contam_pct', ('cluster_ContamPct.tsv',)),
        ('amplitude', ('cluster_Amplitude.tsv',)),
        ('isi_viol', ('cluster_isi_viol.tsv',)),
        ('mergemax', ('cluster_MergeMax.tsv', 'cluster_mergemax.tsv')),
    ):
        for filename in filenames:
            for row in read_rows(directory / filename):
                cid = cluster_id(row)
                if cid in rows:
                    keys = ALIASES.get(feature, (feature,))
                    value = number(row, keys)
                    if not np.isfinite(value):
                        value = number(row, tuple(k for k in row if k not in ('cluster_id', 'id', 'cluster')))
                    if np.isfinite(value):
                        rows[cid][feature] = value
    times = np.load(directory / 'spike_times.npy', mmap_mode='r').reshape(-1)
    if times.size != assignments.size:
        raise ValueError('Spike times and cluster membership lengths differ: ' + str(directory))
    duration, duration_source = recording_duration(directory, parameters(directory), times)
    # Cluster IDs after manual merges/splits need not index the original matrix.
    sim_path = directory / 'similar_templates.npy'
    sim = np.load(sim_path, mmap_mode='r') if sim_path.is_file() else None
    for cid, count in zip(ids, counts):
        row = rows[int(cid)]
        row['n_spikes'] = int(count)
        row['firing_rate'] = float(count) / duration
        if sim is not None and sim.ndim == 2 and sim.shape[0] == sim.shape[1] and 0 <= cid < sim.shape[0]:
            other = [int(x) for x in ids if 0 <= x < sim.shape[0] and x != cid]
            if other:
                values = np.asarray(sim[int(cid), other], dtype=float)
                if np.isfinite(values).any():
                    row['mergemax'] = float(np.nanmax(values))
    matrix = np.array([[rows[int(cid)].get(feature, np.nan) for feature in FEATURES] for cid in ids])
    labels = manual_labels(directory)
    return ids, matrix, labels, {'dataset': str(directory), 'clusters': len(ids),
                                'duration_seconds': duration, 'duration_source': duration_source}


def transform(matrix):
    # Signed log preserves the direction of waveform slopes and velocities.
    return np.sign(matrix) * np.log1p(np.abs(matrix))


def fit_imputer(matrix):
    selected = np.isfinite(matrix).any(axis=0)
    if not selected.any():
        raise ValueError('No finite feature values in training data.')
    selected_matrix = matrix[:, selected]
    medians = np.nanmedian(selected_matrix, axis=0)
    return selected, medians


def impute(matrix, selected, medians):
    values = matrix[:, selected].copy()
    missing = ~np.isfinite(values)
    values[missing] = np.broadcast_to(medians, values.shape)[missing]
    # Missingness is itself recorded, using the same fixed feature schema.
    return np.column_stack((values, missing.astype(float)))


def forest(jobs):
    return RandomForestClassifier(n_estimators=200, max_depth=12, min_samples_leaf=5,
                                  class_weight='balanced_subsample', n_jobs=jobs, random_state=1)


def evaluate(y, scores, threshold):
    predicted = scores >= threshold
    return {'clusters': len(y), 'good': int(np.sum(y)),
            'precision': float(precision_score(y, predicted, zero_division=0)),
            'recall': float(recall_score(y, predicted, zero_division=0)),
            'confusion_matrix_not_good_then_good': confusion_matrix(y, predicted, labels=[0, 1]).tolist(),
            'average_precision': float(average_precision_score(y, scores)),
            'roc_auc': float(roc_auc_score(y, scores)) if len(np.unique(y)) == 2 else None}


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    temporary.replace(path)


def fingerprint(directory):
    return {name: {'size': (directory / name).stat().st_size,
                   'mtime_ns': (directory / name).stat().st_mtime_ns}
            for name in SOURCE_FILES if (directory / name).is_file()}


def train(args):
    if args.model.exists() and not args.overwrite:
        raise ValueError('Model exists; choose another path or pass --overwrite.')
    matrices, targets, groups, sources, weights = [], [], [], [], []
    skipped = []
    roots = list(args.input)
    base_model = getattr(args, 'base_model', None)
    if base_model:
        base = joblib.load(base_model)
        roots.extend(Path(item['dataset']) for item in base['metadata']['datasets'])
    datasets = discover(roots)
    for directory in datasets:
        if not (directory / 'cluster_group.tsv').is_file():
            skipped.append({'dataset': str(directory), 'reason': 'missing manual cluster_group.tsv'})
            continue
        log('TRAIN READ ' + str(directory))
        ids, matrix, labels, info = extract(directory)
        feedback_labels = training_feedback(directory)
        if feedback_labels is not None:
            labels = feedback_labels
        keep = np.array([labels.get(int(cid)) in ('good', 'mua', 'noise') for cid in ids])
        if not keep.any():
            skipped.append({'dataset': str(directory), 'reason': 'no good/mua/noise labels on active clusters'})
            continue
        target = np.array([int(labels.get(int(cid)) == 'good') for cid in ids[keep]])
        matrices.append(matrix[keep])
        targets.append(target)
        weights.extend([getattr(args, 'feedback_weight', 2.0) if feedback_labels is not None else 1.0] * len(target))
        group = recording_group(directory)
        groups.extend([group] * len(target))
        info.update({'recording_group': group, 'labels': dict(Counter(labels.get(int(cid)) for cid in ids[keep])),
                     'label_source': 'saved_explicit_feedback' if feedback_labels is not None else 'legacy_cluster_group',
                     'files': fingerprint(directory)})
        sources.append(info)
    if not matrices:
        raise ValueError('No manually labeled active clusters found.')
    raw = np.vstack(matrices)
    x, y, groups = transform(raw), np.concatenate(targets), np.asarray(groups)
    weights = np.asarray(weights)
    if len(np.unique(y)) != 2:
        raise ValueError('Training needs both good and mua/noise examples.')
    validation = {'status': 'unavailable', 'reason': 'At least two recording groups with both training classes required.'}
    if len(np.unique(groups)) >= 2:
        splitter = GroupShuffleSplit(n_splits=30, test_size=args.validation_fraction, random_state=1)
        for training, heldout in splitter.split(x, y, groups):
            if len(np.unique(y[training])) != 2 or len(np.unique(y[heldout])) != 2:
                continue
            selected, medians = fit_imputer(x[training])
            model = forest(args.jobs).fit(impute(x[training], selected, medians), y[training], sample_weight=weights[training])
            scores = model.predict_proba(impute(x[heldout], selected, medians))[:, 1]
            validation = {'status': 'evaluated', 'method': 'recording-group holdout (all shanks kept together)',
                          'training_recordings': sorted(set(groups[training])),
                          'heldout_recordings': sorted(set(groups[heldout])),
                          'threshold': args.threshold, **evaluate(y[heldout], scores, args.threshold)}
            break
    selected, medians = fit_imputer(x)
    model = forest(args.jobs).fit(impute(x, selected, medians), y, sample_weight=weights)
    metadata = {'format_version': 1, 'program_version': VERSION, 'sklearn_version': sklearn_version,
                'created_utc': datetime.now(timezone.utc).isoformat(), 'features': list(FEATURES),
                'used_features': [f for f, chosen in zip(FEATURES, selected) if chosen],
                'threshold': args.threshold, 'training_clusters': len(y), 'training_good': int(y.sum()),
                'training_other': int(len(y) - y.sum()), 'datasets': sources, 'skipped': skipped,
                'validation': validation,
                'label_policy': 'Saved explicit feedback for tracked datasets; cluster_group.tsv only for legacy curated datasets; no KSLabel fallback.',
                'base_model': str(base_model) if base_model else None,
                'feedback_weight': getattr(args, 'feedback_weight', 2.0),
                'curation_note': 'User-selected training roots must contain manually reviewed curation; file existence alone cannot prove manual review.',
                'feature_importance': dict(zip(
                    [f for f, v in zip(FEATURES, selected) if v] +
                    [f + '_missing' for f, v in zip(FEATURES, selected) if v],
                    model.feature_importances_.tolist()))}
    artifact = {'metadata': metadata, 'model': model, 'selected': selected, 'medians': medians}
    args.model.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.model.with_name(args.model.name + '.tmp')
    joblib.dump(artifact, temporary)
    temporary.replace(args.model)
    write_json(args.model.with_suffix(args.model.suffix + '.json'), metadata)
    log('SAVED ' + str(args.model))
    log(json.dumps({'training_clusters': len(y), 'validation': validation}, ensure_ascii=False))


def write_tsv(path, columns, rows):
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(columns)
        writer.writerows(rows)
    temporary.replace(path)


def predict(args):
    # Only load locally trained/trusted joblib models; pickle-based formats execute code.
    artifact = joblib.load(args.model)
    metadata = artifact['metadata']
    if metadata['format_version'] != 1 or metadata['features'] != list(FEATURES):
        raise ValueError('Unsupported model feature schema.')
    if metadata['sklearn_version'] != sklearn_version:
        raise ValueError('Use the training scikit-learn version: ' + metadata['sklearn_version'])
    threshold = metadata['threshold'] if args.threshold is None else args.threshold
    datasets = discover(args.input)
    filenames = ('prephy_good_predictions.tsv', 'cluster_prephy_good_score.tsv',
                 'cluster_prephy_good_prediction.tsv', 'prephy_good_prediction.json')
    for directory in datasets:
        existing = [name for name in filenames if (directory / name).exists()]
        if existing and not args.overwrite:
            raise ValueError('Prediction files already exist in ' + str(directory) + '; use --overwrite to replace predictions only.')
    model_hash = hashlib.sha256(args.model.read_bytes()).hexdigest()
    summary = []
    for directory in datasets:
        log('PREDICT ' + str(directory))
        ids, raw, labels, info = extract(directory)
        coverage = np.isfinite(raw[:, artifact['selected']]).mean(axis=1)
        probability = artifact['model'].predict_proba(impute(transform(raw), artifact['selected'], artifact['medians']))[:, 1]
        predictions = np.where(coverage < args.min_feature_fraction, 'unknown',
                               np.where(probability >= threshold, 'good', 'not_good'))
        write_tsv(directory / filenames[0],
                  ['cluster_id', 'prephy_good_score', 'prephy_good_prediction', 'prephy_feature_fraction', 'prephy_existing_group'],
                  [(int(cid), '' if label == 'unknown' else format(float(score), '.8f'), str(label), format(float(fraction), '.4f'), labels.get(int(cid), ''))
                   for cid, score, label, fraction in zip(ids, probability, predictions, coverage)])
        write_tsv(directory / filenames[1], ['cluster_id', 'prephy_good_score'],
                  [(int(cid), '' if label == 'unknown' else format(float(score), '.8f'))
                   for cid, score, label in zip(ids, probability, predictions)])
        write_tsv(directory / filenames[2], ['cluster_id', 'prephy_good_prediction'],
                  [(int(cid), str(label)) for cid, label in zip(ids, predictions)])
        info.update({'counts': dict(Counter(str(x) for x in predictions)), 'model_sha256': model_hash,
                     'membership_sha256': array_hash(np.load(directory / 'spike_clusters.npy', mmap_mode='r')),
                     'model': str(args.model.resolve()), 'threshold': threshold,
                     'min_feature_fraction': args.min_feature_fraction, 'source_files': fingerprint(directory),
                     'created_utc': datetime.now(timezone.utc).isoformat(),
                     'training_dataset': str(directory).casefold() in {s['dataset'].casefold() for s in metadata['datasets']}})
        write_json(directory / filenames[3], info)
        summary.append(info)
        log(json.dumps(info['counts']))
    if args.report:
        write_json(args.report, {'datasets': summary})
    log('DONE: ' + str(len(summary)) + ' datasets. Curation labels were not changed.')


def unit_fraction(value):
    value = float(value)
    if not 0 < value < 1:
        raise argparse.ArgumentTypeError('Expected a value strictly between 0 and 1.')
    return value


def main():
    # Use a consistent encoding for redirected logs on Japanese Windows.
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, 'reconfigure'):
            stream.reconfigure(encoding='utf-8', errors='backslashreplace')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=VERSION)
    commands = parser.add_subparsers(dest='command', required=True)
    learning = commands.add_parser('train', help='Train using user-selected manually curated datasets.')
    prediction = commands.add_parser('predict', help='Write prediction TSVs without changing cluster_group.tsv.')
    for command in (learning, prediction):
        command.add_argument('--input', type=Path, action='append', required=True,
                             help='Dataset or parent directory; repeat to include multiple roots.')
        command.add_argument('--model', type=Path, required=True)
        command.add_argument('--overwrite', action='store_true')
    learning.add_argument('--threshold', type=unit_fraction, default=0.5)
    learning.add_argument('--validation-fraction', type=unit_fraction, default=0.25)
    learning.add_argument('--jobs', type=int, default=1)
    learning.add_argument('--base-model', type=Path, help='Retrain from this model\'s source datasets plus --input feedback directories.')
    learning.add_argument('--feedback-weight', type=float, default=2.0, help='Weight for explicitly reviewed feedback examples (default 2).')
    learning.set_defaults(function=train)
    prediction.add_argument('--threshold', type=unit_fraction)
    prediction.add_argument('--min-feature-fraction', type=unit_fraction, default=0.5)
    prediction.add_argument('--report', type=Path)
    prediction.set_defaults(function=predict)
    args = parser.parse_args()
    if args.model.suffix.lower() != '.joblib':
        parser.error('--model must have the .joblib extension.')
    if getattr(args, 'report', None) and args.report.suffix.lower() != '.json':
        parser.error('--report must have the .json extension.')
    if getattr(args, 'jobs', 1) == 0:
        parser.error('--jobs cannot be zero.')
    if not np.isfinite(getattr(args, 'feedback_weight', 2.0)) or getattr(args, 'feedback_weight', 2.0) <= 0:
        parser.error('--feedback-weight must be positive and finite.')
    try:
        args.function(args)
    except (ValueError, OSError) as error:
        parser.exit(1, 'ERROR: ' + str(error) + '\n')


if __name__ == '__main__':
    main()
