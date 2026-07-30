r"""Build a Rastermap matrix from phy-curated good Kilosort clusters.

Example:
    python run_good_clusters_rastermap.py . --out ./rastermap \
        --bin-ms 100 --run-rastermap
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class Dataset:
    path: Path
    sample_rate: float
    good_cluster_ids: np.ndarray
    spike_times_seconds: np.ndarray
    spike_clusters: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a good-cluster activity matrix from imec*_ks4 folders "
            "and optionally run Rastermap."
        )
    )
    parser.add_argument(
        "roots", nargs="+", type=Path,
        help="Root directories to scan recursively for imec*_ks4 datasets.",
    )
    parser.add_argument(
        "--out", type=Path, default=Path("rastermap_good_clusters"),
        help="Output directory.",
    )
    parser.add_argument(
        "--dataset-glob", default="imec*_ks4",
        help="Dataset directory name pattern. Default: imec*_ks4.",
    )
    parser.add_argument("--bin-ms", type=float, default=100.0)
    parser.add_argument(
        "--sample-rate", type=float,
        help="Override the sample rate when params.py cannot be read.",
    )
    parser.add_argument("--min-spikes", type=int, default=1)
    parser.add_argument(
        "--normalize", choices=("zscore", "sqrt_zscore", "none"),
        default="sqrt_zscore",
    )
    parser.add_argument(
        "--max-neurons", type=int, default=0,
        help="Keep the most active N good clusters; 0 keeps all.",
    )
    parser.add_argument(
        "--max-time-bins", type=int, default=0,
        help="Uniformly sample N time bins; 0 keeps all.",
    )
    parser.add_argument("--run-rastermap", action="store_true")
    parser.add_argument("--n-pcs", type=int, default=200)
    parser.add_argument("--n-clusters", type=int, default=100)
    parser.add_argument("--locality", type=float, default=0.75)
    parser.add_argument("--time-lag-window", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=0)
    return parser.parse_args()


def find_datasets(roots: Iterable[Path], pattern: str) -> list[Path]:
    datasets: list[Path] = []
    for root in roots:
        root = root.expanduser()
        if not root.exists():
            print(f"missing root: {root}", file=sys.stderr)
            continue
        candidates = [root] if root.match(pattern) else root.rglob(pattern)
        for path in candidates:
            if not path.is_dir() or "_orig" in path.name.lower():
                continue
            required = ("spike_times.npy", "spike_clusters.npy")
            if all((path / name).exists() for name in required):
                datasets.append(path.resolve())
    return sorted(set(datasets))


def read_cluster_groups(dataset_dir: Path) -> dict[int, str]:
    path = dataset_dir / "cluster_group.tsv"
    if not path.exists():
        return {}
    groups: dict[int, str] = {}
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cluster = row.get("cluster_id") or row.get("id") or row.get("cluster")
            group = row.get("group") or row.get("KSLabel") or row.get("label")
            if cluster in (None, "") or group in (None, ""):
                continue
            try:
                groups[int(cluster)] = str(group).strip().lower()
            except ValueError:
                continue
    return groups


def read_sample_rate(dataset_dir: Path, override: float | None) -> float:
    if override is not None:
        if override <= 0:
            raise ValueError("--sample-rate must be positive")
        return float(override)
    params = dataset_dir / "params.py"
    if not params.exists():
        return 30000.0
    for line in params.read_text(encoding="utf-8", errors="ignore").splitlines():
        if "=" not in line:
            continue
        key, raw = (part.strip() for part in line.split("=", 1))
        if key not in {"sample_rate", "fs", "sampling_rate"}:
            continue
        try:
            return float(ast.literal_eval(raw))
        except (ValueError, SyntaxError):
            try:
                return float(raw)
            except ValueError:
                pass
    return 30000.0


def load_dataset(path: Path, sample_rate: float | None) -> Dataset | None:
    groups = read_cluster_groups(path)
    good_ids = sorted(cluster for cluster, label in groups.items() if label == "good")
    if not good_ids:
        return None
    spike_times_raw = np.asarray(np.load(path / "spike_times.npy")).reshape(-1)
    spike_clusters = (
        np.asarray(np.load(path / "spike_clusters.npy"))
        .reshape(-1).astype(np.int64, copy=False)
    )
    if spike_times_raw.size != spike_clusters.size:
        raise ValueError(f"{path}: spike_times and spike_clusters have different lengths")
    rate = read_sample_rate(path, sample_rate)
    spike_times_float = spike_times_raw.astype(np.float64, copy=False)
    if np.issubdtype(spike_times_raw.dtype, np.integer):
        spike_times_seconds = spike_times_float / rate
    elif spike_times_float.size and np.nanmax(spike_times_float) > 1e6:
        spike_times_seconds = spike_times_float / rate
    else:
        spike_times_seconds = spike_times_float
    return Dataset(
        path=path,
        sample_rate=rate,
        good_cluster_ids=np.asarray(good_ids, dtype=np.int64),
        spike_times_seconds=spike_times_seconds,
        spike_clusters=spike_clusters,
    )


def build_activity_matrix(
    datasets: list[Dataset],
    bin_ms: float,
    min_spikes: int,
    max_neurons: int,
    max_time_bins: int,
) -> tuple[np.ndarray, list[dict[str, object]], np.ndarray]:
    bin_seconds = bin_ms / 1000.0
    if bin_seconds <= 0 or min_spikes < 1:
        raise ValueError("--bin-ms and --min-spikes must be positive")
    max_time = max(
        (
            float(np.nanmax(dataset.spike_times_seconds))
            for dataset in datasets if dataset.spike_times_seconds.size
        ),
        default=0.0,
    )
    n_bins = max(1, int(math.ceil(max_time / bin_seconds)) + 1)
    rows: list[np.ndarray] = []
    metadata: list[dict[str, object]] = []
    for dataset in datasets:
        for cluster_id in dataset.good_cluster_ids:
            mask = dataset.spike_clusters == int(cluster_id)
            n_spikes = int(mask.sum())
            if n_spikes < min_spikes:
                continue
            bins = np.floor(
                dataset.spike_times_seconds[mask] / bin_seconds
            ).astype(np.int64, copy=False)
            bins = bins[(bins >= 0) & (bins < n_bins)]
            rows.append(np.bincount(bins, minlength=n_bins).astype(np.float32))
            metadata.append({
                "dataset": str(dataset.path),
                "dataset_name": dataset.path.name,
                "cluster_id": int(cluster_id),
                "n_spikes": n_spikes,
                "sample_rate": dataset.sample_rate,
            })
    if not rows:
        raise RuntimeError("No good clusters passed the filters.")
    matrix = np.vstack(rows).astype(np.float32, copy=False)
    if max_neurons and matrix.shape[0] > max_neurons:
        order = np.argsort(matrix.sum(axis=1))[::-1][:max_neurons]
        matrix = matrix[order]
        metadata = [metadata[int(index)] for index in order]
    time_seconds = np.arange(matrix.shape[1], dtype=np.float32) * bin_seconds
    if max_time_bins and matrix.shape[1] > max_time_bins:
        indices = np.linspace(
            0, matrix.shape[1] - 1, max_time_bins
        ).astype(np.int64)
        matrix = matrix[:, indices]
        time_seconds = time_seconds[indices]
    return matrix, metadata, time_seconds


def normalize_activity(matrix: np.ndarray, mode: str) -> np.ndarray:
    normalized = matrix.astype(np.float32, copy=True)
    if mode == "none":
        return normalized
    if mode == "sqrt_zscore":
        normalized = np.sqrt(normalized)
    standard_deviation = normalized.std(axis=1, keepdims=True)
    standard_deviation[standard_deviation < 1e-6] = 1.0
    return (
        normalized - normalized.mean(axis=1, keepdims=True)
    ) / standard_deviation


def fit_rastermap(matrix: np.ndarray, args: argparse.Namespace) -> dict[str, np.ndarray]:
    try:
        from rastermap import Rastermap
    except ImportError as error:
        raise RuntimeError(
            "rastermap is not installed; run: python -m pip install rastermap"
        ) from error
    if min(matrix.shape) < 2:
        raise RuntimeError("Rastermap needs at least two clusters and two time bins.")
    kwargs = {
        "n_PCs": min(args.n_pcs, max(1, min(matrix.shape) - 1)),
        "n_clusters": min(args.n_clusters, max(2, matrix.shape[0])),
        "locality": args.locality,
        "time_lag_window": args.time_lag_window,
        "random_state": args.random_state,
    }
    try:
        model = Rastermap(**kwargs).fit(matrix)
    except TypeError:
        kwargs.pop("random_state")
        model = Rastermap(**kwargs).fit(matrix)
    result = {
        "isort": np.asarray(model.isort),
        "embedding": np.asarray(model.embedding),
    }
    if hasattr(model, "X_embedding"):
        result["X_embedding"] = np.asarray(model.X_embedding)
    return result


def save_preview(out: Path, matrix: np.ndarray, title: str, name: str) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; preview was skipped", file=sys.stderr)
        return
    figure, axis = plt.subplots(figsize=(12, 6))
    axis.imshow(
        matrix, aspect="auto", cmap="gray_r",
        vmin=np.percentile(matrix, 1.0), vmax=np.percentile(matrix, 99.5),
    )
    axis.set(xlabel="time bin", ylabel="good cluster", title=title)
    figure.tight_layout()
    figure.savefig(out / name, dpi=160)
    plt.close(figure)


def main() -> int:
    args = parse_args()
    out = args.out.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    paths = find_datasets(args.roots, args.dataset_glob)
    print(f"found {len(paths)} candidate datasets")
    datasets = []
    for path in paths:
        dataset = load_dataset(path, args.sample_rate)
        if dataset is not None:
            datasets.append(dataset)
            print(f"loaded {path}: good_clusters={len(dataset.good_cluster_ids)}")
    if not datasets:
        raise RuntimeError("No datasets with good clusters were found.")
    matrix, metadata, time_seconds = build_activity_matrix(
        datasets, args.bin_ms, args.min_spikes,
        args.max_neurons, args.max_time_bins,
    )
    normalized = normalize_activity(matrix, args.normalize)
    np.save(out / "spks_good.npy", matrix)
    np.save(out / "spks_good_normalized.npy", normalized)
    np.save(out / "time_seconds.npy", time_seconds)
    (out / "cluster_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    summary = {
        "n_datasets": len(datasets),
        "n_good_clusters": len(metadata),
        "n_time_bins": int(matrix.shape[1]),
        "bin_ms": args.bin_ms,
        "normalize": args.normalize,
        "roots": [str(path) for path in args.roots],
    }
    (out / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    save_preview(
        out, normalized, "Good clusters before Rastermap sorting",
        "good_clusters_unsorted.png",
    )
    if args.run_rastermap:
        result = fit_rastermap(normalized, args)
        np.savez(out / "rastermap_result.npz", **result)
        order = result["isort"].astype(np.int64, copy=False)
        np.save(out / "rastermap_isort.npy", order)
        sorted_metadata = [metadata[int(index)] for index in order]
        (out / "cluster_metadata_rastermap_sorted.json").write_text(
            json.dumps(sorted_metadata, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        display = result.get("X_embedding", normalized[order])
        save_preview(
            out, display, "Good clusters sorted by Rastermap",
            "good_clusters_rastermap.png",
        )
        print(f"rastermap complete: {out / 'rastermap_result.npz'}")
    else:
        print("Rastermap skipped; add --run-rastermap to fit it.")
    print(f"wrote outputs to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
