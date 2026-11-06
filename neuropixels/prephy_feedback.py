"""Saved human feedback and one-time initialization; no model/GUI dependency."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np


def now():
    return datetime.now(timezone.utc).isoformat()


def array_hash(array):
    array = np.asarray(array).reshape(-1).astype('<i8')
    return hashlib.sha256(array.tobytes()).hexdigest()


def signatures(array):
    array = np.asarray(array).reshape(-1)
    order = np.argsort(array, kind='stable')
    ids, counts = np.unique(array, return_counts=True)
    out, start = {}, 0
    for cid, count in zip(ids, counts):
        out[str(int(cid))] = hashlib.sha256(order[start:start + count].astype('<i8').tobytes()).hexdigest()
        start += count
    return out


def read_groups(directory):
    path = Path(directory) / 'cluster_group.tsv'
    if not path.is_file():
        return {}
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return {int(row['cluster_id']): row.get('group', '').strip().lower()
                for row in csv.DictReader(handle, delimiter='\t')}


def write_json(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')
    temporary.replace(path)


class FeedbackStore:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.folder = self.directory / 'prephy_feedback'
        self.path = self.folder / 'state.json'

    def load(self):
        return json.loads(self.path.read_text(encoding='utf-8')) if self.path.is_file() else None

    def initialize(self, assignments):
        """Return assignments to apply once, preserving any existing curation."""
        state = self.load()
        if state:
            if state.get('saved_membership_sha256') or state.get('initialization_policy') != 'prediction_good_else_noise':
                return {}
            # If phy was closed without saving the initial labels, retry on the
            # next opening. A saved/external manual group always takes priority.
            existing = read_groups(self.directory)
            if any(label in ('good', 'mua', 'noise') for label in existing.values()):
                return {}
            if array_hash(assignments) != state.get('initial_membership_sha256'):
                raise ValueError('Membership changed since initial opening; initialization was not repeated.')
            return {int(cid): label for cid, label in state.get('initial_assignments', {}).items()}
        prediction = self.directory / 'prephy_good_predictions.tsv'
        if not prediction.is_file():
            return {}
        active = set(int(x) for x in np.unique(assignments))
        with prediction.open(encoding='utf-8-sig', newline='') as handle:
            predictions = {int(row['cluster_id']): row for row in csv.DictReader(handle, delimiter='\t')}
        existing = read_groups(self.directory)
        curated = any(existing.get(cid) in ('good', 'mua', 'noise') for cid in active)
        policy = 'existing_curation_preserved' if curated else 'prediction_good_else_noise'
        report_path = self.directory / 'prephy_good_prediction.json'
        report = json.loads(report_path.read_text(encoding='utf-8')) if report_path.is_file() else {}
        if not curated:
            if active != set(predictions):
                raise ValueError('Prediction cluster IDs are stale. Run predict again before initialization.')
            expected = report.get('membership_sha256')
            if expected:
                valid = expected == array_hash(assignments)
            else:
                # Compatibility with predictions produced by prephy 1.0.
                stamp = report.get('source_files', {}).get('spike_clusters.npy', {})
                stat = (self.directory / 'spike_clusters.npy').stat()
                valid = stamp.get('mtime_ns') == stat.st_mtime_ns and stamp.get('size') == stat.st_size
            if not valid:
                raise ValueError('Prediction membership provenance is missing or stale. Run predict again.')
        plan = {} if curated else {cid: ('good' if predictions[cid]['prephy_good_prediction'] == 'good' else 'noise') for cid in active}
        self.folder.mkdir(exist_ok=True)
        group_file = self.directory / 'cluster_group.tsv'
        if group_file.is_file():
            shutil.copyfile(group_file, self.folder / 'cluster_group.before_initialization.tsv')
        shutil.copyfile(prediction, self.folder / 'predictions_initial.tsv')
        if report_path.is_file():
            shutil.copyfile(report_path, self.folder / 'prediction_initial.json')
        np.save(self.folder / 'membership_initial.npy', np.asarray(assignments).reshape(-1))
        state = {'version': 1, 'created': now(), 'initialization_applied': not bool(plan),
                 'initialization_policy': policy, 'initial_assignments': {str(k): v for k, v in plan.items()},
                 'initial_membership_sha256': array_hash(assignments),
                 'reviewed': {}, 'complete': False, 'save_number': 0,
                 'saved_membership_sha256': None}
        write_json(self.path, state)
        return plan

    def initialized(self):
        state = self.load()
        if state:
            state['initialization_applied'] = True
            write_json(self.path, state)

    def commit(self, assignments, groups, reviewed_ids=(), complete=False, actions=()):
        state = self.load()
        if state is None:
            return
        assignments = np.asarray(assignments).reshape(-1)
        # Called after phy's own save callback. Only persist feedback for a
        # successfully saved membership and its matching group file.
        saved = np.load(self.directory / 'spike_clusters.npy', mmap_mode='r').reshape(-1)
        digest = array_hash(assignments)
        if array_hash(saved) != digest:
            raise ValueError('phy has not saved this spike membership; feedback was not committed.')
        active = set(int(x) for x in np.unique(assignments))
        saved_groups = read_groups(self.directory)
        groups = {int(k): str(v).lower() for k, v in groups.items()}
        if any(saved_groups.get(cid, 'unsorted') != groups.get(cid, 'unsorted') for cid in active):
            raise ValueError('phy has not saved these group labels; feedback was not committed.')
        sigs = signatures(assignments)
        previous = state.get('reviewed', {})
        reviewed = {cid: item for cid, item in previous.items()
                    if cid in sigs and item['membership_signature'] == sigs[cid]
                    and item['label'] == groups.get(int(cid))}
        explicit = active if complete else active.intersection(int(x) for x in reviewed_ids)
        if complete and any(groups.get(cid) not in ('good', 'mua', 'noise') for cid in active):
            raise ValueError('Assign good/mua/noise to all clusters before confirming completion.')
        for cid in explicit:
            if groups.get(cid) in ('good', 'mua', 'noise'):
                reviewed[str(cid)] = {'label': groups[cid], 'membership_signature': sigs[str(cid)], 'saved_utc': now()}
        number = int(state.get('save_number', 0)) + 1
        stem = 'save_' + str(number).zfill(6)
        membership_file = state.get('membership_snapshot', 'membership_initial.npy')
        if state.get('saved_membership_sha256') != digest:
            membership_file = stem + '_spike_clusters.npy'
            np.save(self.folder / membership_file, assignments)
        feature_file = self.directory / 'cluster_info.tsv'
        if feature_file.is_file():
            shutil.copyfile(feature_file, self.folder / (stem + '_cluster_info.tsv'))
        event = {'time': now(), 'save_number': number, 'actions': list(actions),
                 'explicitly_reviewed_ids': sorted(explicit), 'confirmed_complete': bool(complete),
                 'groups': {str(k): v for k, v in groups.items()},
                 'prediction_snapshot': 'predictions_initial.tsv',
                 'features_snapshot': stem + '_cluster_info.tsv' if feature_file.is_file() else None,
                 'membership_snapshot': membership_file, 'membership_sha256': digest}
        with (self.folder / 'history.jsonl').open('a', encoding='utf-8') as handle:
            handle.write(json.dumps(event, ensure_ascii=False) + '\n')
        state.update({'reviewed': reviewed, 'complete': bool(active) and len(reviewed) == len(active),
                      'save_number': number, 'saved_membership_sha256': digest,
                      'membership_snapshot': membership_file, 'updated': now()})
        write_json(self.path, state)


def training_feedback(directory):
    """None for legacy curation; otherwise only saved, explicitly reviewed labels."""
    store = FeedbackStore(directory)
    state = store.load()
    if state is None:
        return None
    expected = state.get('saved_membership_sha256')
    current = np.load(Path(directory) / 'spike_clusters.npy', mmap_mode='r')
    if not expected or array_hash(current) != expected:
        return {}
    labels = read_groups(directory)
    return {int(cid): item['label'] for cid, item in state.get('reviewed', {}).items()
            if labels.get(int(cid)) == item['label']}
