from __future__ import annotations

from dataclasses import dataclass, asdict
from itertools import product
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Any

import numpy as np
import pandas as pd

from sklearn.decomposition import PCA as SklearnPCA
from deeptime.decomposition import TICA, VAMP
from deeptime.clustering import KMeans
from deeptime.markov import TransitionCountEstimator
from deeptime.markov.msm import MaximumLikelihoodMSM


Array = np.ndarray
Trajs = List[Array]


def concatenate(trajs: Trajs) -> Array:
    """
    Concatenates a list of trajectories into a single array.
    """
    if len(trajs) == 0:
        raise ValueError("Cannot concatenate an empty list of trajectories.")
    return np.concatenate(trajs, axis=0)


def split_like_concatenated(concat: Array, reference_trajs: Trajs) -> Trajs:
    """
    Splits a concatenated array back into a list of trajectories with the same
    lengths as the reference trajectories.
    """
    lengths = [len(x) for x in reference_trajs]
    out = []
    start = 0
    for n in lengths:
        out.append(concat[start:start + n])
        start += n
    return out


@dataclass(frozen=True)
class PipelineParams:
    """
    Parameters defining one MSM construction pipeline.

    frozen=True makes the dataclass immutable after creation. This is useful for
    grid searches because each parameter set is treated as a fixed, hashable-like
    configuration rather than something accidentally modified in place.
    """
    dimred_method: str = "tica"
    dimred_lagtime: int = 10
    dim: int = 5

    n_clusters: int = 100
    clustering_method: str = "kmeans"

    msm_lagtime: int = 10
    reversible: bool = True
    count_mode: str = "sliding-effective"

    vamp_r: Any = 2
    vamp_cv_splits: int = 5


class DimReducerFactory:
    """
    Registry-style dimensionality reduction factory.
    Add new methods by registering a new constructor.
    """

    def __init__(self):
        self._registry: Dict[str, Callable[[PipelineParams], Any]] = {}

    def register(self, name: str, constructor: Callable[[PipelineParams], Any]) -> None:
        self._registry[name] = constructor

    def make(self, params: PipelineParams):
        if params.dimred_method not in self._registry:
            raise ValueError(f"Unknown dimred method: {params.dimred_method}")
        return self._registry[params.dimred_method](params)


dimred_factory = DimReducerFactory()

dimred_factory.register(
    "tica",
    lambda p: TICA(
        lagtime=p.dimred_lagtime,
        dim=p.dim,
    ),
)


class PCAWrapper:
    """
    Small wrapper that gives sklearn PCA a deeptime-like API:
    fit(...).fetch_model().transform(...).
    """

    def __init__(self, dim: int):
        self.model = SklearnPCA(n_components=dim)

    def fit(self, trajs: Trajs):
        X = concatenate(trajs)
        self.model.fit(X)
        return self

    def fetch_model(self):
        return self

    def transform(self, x: Array) -> Array:
        return self.model.transform(x)


dimred_factory.register(
    "pca",
    lambda p: PCAWrapper(dim=p.dim),
)


# Example for registering an additional method:
#
# from deeptime.decomposition import VAMP
# dimred_factory.register(
#     "vamp",
#     lambda p: VAMP(
#         lagtime=p.dimred_lagtime,
#         dim=p.dim,
#     ),
# )


def make_clusterer(params: PipelineParams):
    if params.clustering_method == "kmeans":
        return KMeans(
            n_clusters=params.n_clusters,
            max_iter=50,
            fixed_seed=13,
        )
    raise ValueError(f"Unknown clustering method: {params.clustering_method}")


def fit_transform_dimred(trajs: Trajs, params: PipelineParams) -> Tuple[Trajs, Any]:
    """
    1. Fits the dimensionality reduction model.
    2. Transforms the trajectories.
    3. Returns reduced trajectories and fitted model.
    """
    estimator = dimred_factory.make(params)
    model = estimator.fit(trajs).fetch_model()
    reduced = [model.transform(x) for x in trajs]
    return reduced, model


def fit_transform_clustering(reduced_trajs: Trajs, params: PipelineParams) -> Tuple[List[np.ndarray], Any]:
    """
    1. Fits clustering on concatenated reduced trajectories.
    2. Transforms each reduced trajectory into a discrete trajectory.
    3. Returns dtrajs and fitted clustering model.
    """
    clusterer = make_clusterer(params)

    X = concatenate(reduced_trajs)
    cluster_model = clusterer.fit(X).fetch_model()

    dtrajs = [cluster_model.transform(x).astype(int) for x in reduced_trajs]
    return dtrajs, cluster_model


def fit_msm(dtrajs: List[np.ndarray], params: PipelineParams):
    """
    Fits an MSM model on the discrete trajectories.
    """
    counts_estimator = TransitionCountEstimator(
        lagtime=params.msm_lagtime,
        count_mode=params.count_mode,
    )
    counts = counts_estimator.fit(dtrajs).fetch_model()

    msm_estimator = MaximumLikelihoodMSM(
        reversible=params.reversible,
    )
    msm = msm_estimator.fit(counts).fetch_model()
    return msm


def fit_pipeline(trajs: Trajs, params: PipelineParams) -> Dict[str, Any]:
    """
    Fits the full final pipeline:

        dimred -> clustering -> MLE MSM

    This is used after hyperparameter selection, when fitting the final model on
    the whole dataset.
    """
    reduced_trajs, dimred_model = fit_transform_dimred(trajs, params)
    dtrajs, cluster_model = fit_transform_clustering(reduced_trajs, params)
    msm = fit_msm(dtrajs, params)

    return {
        "params": params,
        "dimred_model": dimred_model,
        "reduced_trajs": reduced_trajs,
        "cluster_model": cluster_model,
        "dtrajs": dtrajs,
        "msm": msm,
    }


# -----------------------------------------------------------------------------
# Cross-validated variational scoring
# -----------------------------------------------------------------------------


def one_hot_dtraj(dtraj: np.ndarray, n_states: int) -> np.ndarray:
    """
    Converts a discrete trajectory into one-hot MSM indicator functions.
    """
    dtraj = np.asarray(dtraj, dtype=int)

    if len(dtraj) == 0:
        return np.zeros((0, n_states), dtype=np.float64)

    if np.any(dtraj < 0) or np.any(dtraj >= n_states):
        raise ValueError("Discrete trajectory contains state labels outside [0, n_states).")

    X = np.zeros((len(dtraj), n_states), dtype=np.float64)
    X[np.arange(len(dtraj)), dtraj] = 1.0
    return X


def one_hot_dtrajs(dtrajs: List[np.ndarray], n_states: int) -> List[np.ndarray]:
    return [one_hot_dtraj(dtraj, n_states) for dtraj in dtrajs]


def split_trajs_into_blocks(trajs: Trajs, n_splits: int) -> List[Trajs]:
    """
    Splits each trajectory into contiguous blocks.

    Fold k contains block k from each original trajectory. This avoids random
    frame shuffling and is more appropriate for time-correlated MD data.
    """
    if n_splits < 2:
        raise ValueError("n_splits must be at least 2.")

    folds: List[Trajs] = [[] for _ in range(n_splits)]

    for traj in trajs:
        if len(traj) < n_splits:
            raise ValueError(
                f"Trajectory of length {len(traj)} is shorter than n_splits={n_splits}."
            )

        edges = np.linspace(0, len(traj), n_splits + 1, dtype=int)

        for k in range(n_splits):
            block = traj[edges[k]:edges[k + 1]]
            if len(block) > 0:
                folds[k].append(block)

    return folds


def fit_dimred_and_cluster(train_trajs: Trajs, params: PipelineParams) -> Dict[str, Any]:
    """
    Fits only dimred + clustering on training data.

    The final MSM is not fitted here because the VAMP/VAC-style CV score should
    evaluate the discrete-state basis, not the training-set MSM eigenvalues.
    """
    reduced_train, dimred_model = fit_transform_dimred(train_trajs, params)
    dtrajs_train, cluster_model = fit_transform_clustering(reduced_train, params)

    return {
        "dimred_model": dimred_model,
        "reduced_train": reduced_train,
        "cluster_model": cluster_model,
        "dtrajs_train": dtrajs_train,
    }


def transform_with_fitted_pipeline(
    trajs: Trajs,
    dimred_model: Any,
    cluster_model: Any,
) -> List[np.ndarray]:
    """
    Applies already-fitted dimred and clustering models to new trajectories.
    """
    reduced = [dimred_model.transform(x) for x in trajs]
    dtrajs = [cluster_model.transform(x).astype(int) for x in reduced]
    return dtrajs


def score_params_vamp_cv(
    trajs: Trajs,
    params: PipelineParams,
    random_state: Optional[int] = 42,
    n_jobs: int = 1,
) -> Dict[str, Any]:
    """
    Cross-validated variational score for the discrete MSM basis.

    Logic:
    1. Split trajectories into contiguous folds.
    2. Fit dimred + clustering on the training folds only.
    3. Assign both train and test folds to the same learned clusters.
    4. One-hot encode the discrete trajectories.
    5. Fit Deeptime VAMP on the training one-hot trajectories.
    6. Score the learned slow subspace on the test one-hot trajectories.

    For equilibrium reversible dynamics, this VAMP-2 score may be interpreted as
    the VAC-2 variational score of the MSM indicator basis.

    Parameters random_state and n_jobs are kept for API compatibility with your
    previous version. The current deterministic blocked CV implementation does
    not use them, except that KMeans has a fixed seed in make_clusterer().
    """
    _ = random_state
    _ = n_jobs

    folds = split_trajs_into_blocks(trajs, params.vamp_cv_splits)
    scores = []

    for k in range(params.vamp_cv_splits):
        test_trajs = folds[k]
        train_trajs = [
            traj
            for j, fold in enumerate(folds)
            if j != k
            for traj in fold
        ]

        # Fit transformations/state definition on train only.
        fitted = fit_dimred_and_cluster(train_trajs, params)

        dtrajs_train = fitted["dtrajs_train"]
        dtrajs_test = transform_with_fitted_pipeline(
            test_trajs,
            fitted["dimred_model"],
            fitted["cluster_model"],
        )

        # Convert discrete cluster labels into MSM indicator functions.
        Z_train = one_hot_dtrajs(dtrajs_train, params.n_clusters)
        Z_test = one_hot_dtrajs(dtrajs_test, params.n_clusters)

        vamp = VAMP(lagtime=params.msm_lagtime)
        vamp_model = vamp.fit(Z_train).fetch_model()

        score = vamp_model.score(Z_test, r=params.vamp_r)
        scores.append(score)

    scores = np.asarray(scores, dtype=float)

    return {
        **asdict(params),
        "vamp_mean": float(np.mean(scores)),
        "vamp_std": float(np.std(scores)),
        "vamp_scores": scores,
    }


def expand_grid(grid: Dict[str, Sequence[Any]]) -> Iterable[PipelineParams]:
    """
    Expands a grid of parameters into PipelineParams instances.
    """
    keys = list(grid.keys())
    for values in product(*[grid[k] for k in keys]):
        kwargs = dict(zip(keys, values))
        yield PipelineParams(**kwargs)


def _coerce_best_params_from_row(row: pd.Series) -> PipelineParams:
    """
    Converts a DataFrame row back into PipelineParams while avoiding common
    pandas dtype issues, such as np.int64 instead of int.
    """
    kwargs = {}
    for field_name, field_def in PipelineParams.__dataclass_fields__.items():
        value = row[field_name]

        if pd.isna(value) and field_def.default is None:
            value = None

        if field_name in {
            "dimred_lagtime",
            "dim",
            "n_clusters",
            "msm_lagtime",
            "vamp_cv_splits",
        }:
            value = int(value)

        elif field_name == "reversible":
            value = bool(value)

        kwargs[field_name] = value

    return PipelineParams(**kwargs)


def grid_search(
    trajs: Trajs,
    grid: Dict[str, Sequence[Any]],
    random_state: Optional[int] = 42,
    n_jobs: int = 1,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """
    Scores every parameter combination and fits the best pipeline on all data.

    Returns
    -------
    results : pd.DataFrame
        Sorted table of scores.

    best_model : dict
        Output of fit_pipeline(trajs, best_params).
    """
    rows = []

    for params in expand_grid(grid):
        print(f"Scoring: {params}")
        try:
            row = score_params_vamp_cv(
                trajs=trajs,
                params=params,
                random_state=random_state,
                n_jobs=n_jobs,
            )
            rows.append(row)

        except Exception as exc:
            failed = asdict(params)
            failed.update(
                {
                    "vamp_mean": np.nan,
                    "vamp_std": np.nan,
                    "error": repr(exc),
                }
            )
            rows.append(failed)

    results = pd.DataFrame(rows)
    results = results.sort_values("vamp_mean", ascending=False, na_position="last")

    #if results["vamp_mean"].notna().sum() == 0:
    #    raise RuntimeError(
    #        "All parameter combinations failed. Inspect the 'error' column of the results DataFrame."
    #    )

    if results["vamp_mean"].notna().sum() != 0:
        best_row = results[results["vamp_mean"].notna()].iloc[0]
        best_params = _coerce_best_params_from_row(best_row)
        best_model = fit_pipeline(trajs, best_params)

        return results, best_model
    else:
        print("Warning: All parameter combinations failed. Returning NaN scores and error messages.")
        return results, None
