#!/usr/bin/env python3
"""Independently reconstruct capitulum-only Azami morphometric coupling.

Inputs: pinned numerical archives azami Actions 9612943217 and 9633419268.
NOT new pollination, mechanical-spine, selection, genetic, or coevolution data.
The script does not modify Azami's frozen GEB manuscript.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

INPUTS = {
    "traits_zip": "101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e",
    "environment_zip": "d7c0c466f55b67695d06ae46c21a6452dbe6cfd92a52db8042caa200429e97f4",
}
C = {
    "presentation_angle": ("scalar", ["orientation_image_vertical_angle"], [1]),
    "floral_lightness": ("scalar", ["corolla_lab_lightness"], [1]),
    "floral_chroma": ("scalar", ["corolla_lab_chroma"], [1]),
    "floral_hue": ("joint", ["corolla_hue_sin", "corolla_hue_cos"], []),
    "head_elongation": ("scalar", ["capitulum_outline_aspect_ratio"], [1]),
    "head_compactness": ("scalar", [
        "capitulum_outline_circularity", "capitulum_outline_solidity",
        "capitulum_width_profile_cv"], [1, 1, -1]),
    "involucre_form": ("joint", [
        "involucre_length_width_ratio", "involucre_apical_taper_ratio",
        "involucre_basal_taper_ratio"], []),
    "projection_prominence": ("scalar", [
        "bract_projection_roughness", "bract_projection_p95",
        "bract_projection_maximum", "bract_spread_fraction"], [1, 1, 1, 1]),
    "projection_pattern": ("joint", [
        "bract_projection_peak_density", "bract_projection_asymmetry"], []),
}
PRIMARY = [
    ("presentation_angle", "head_elongation"),
    ("floral_chroma", "projection_prominence"),
    ("presentation_angle", "projection_prominence"),
    ("projection_prominence", "projection_pattern"),
    ("presentation_angle", "involucre_form"),
    ("floral_lightness", "floral_hue"),
]
QC = ["latitude", "chelsa_bio12", "chelsa_rsds_mean",
      "min_dimension", "sharpness", "mask_quality"]


def verify_zip(path: Path, key: str) -> None:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for piece in iter(lambda: f.read(1 << 20), b""):
            h.update(piece)
    if h.hexdigest() != INPUTS[key]:
        raise ValueError(f"Invalid pinned source ZIP hash: {key}")


def rv(x, y, weights=None):
    if weights is None:
        weights = np.ones(len(x))
    w = np.asarray(weights, dtype=float)
    w /= w.sum()
    cx = x - (x * w[:, None]).sum(0)
    cy = y - (y * w[:, None]).sum(0)
    xx = (cx * w[:, None]).T @ cx
    yy = (cy * w[:, None]).T @ cy
    xy = (cx * w[:, None]).T @ cy
    den = np.sqrt(np.sum(xx ** 2) * np.sum(yy ** 2))
    return float(np.sum(xy ** 2) / den) if den > 1e-15 else np.nan


def batch_rv(xx, yy, xy):
    return np.sum(xy ** 2, axis=(-2, -1)) / np.sqrt(
        np.sum(xx ** 2, axis=(-2, -1)) * np.sum(yy ** 2, axis=(-2, -1)))


def partial_corr(x, y, matrix):
    m = (matrix - matrix.mean(0)) / np.where(
        matrix.std(0) > 0, matrix.std(0), 1)
    z = np.column_stack([np.ones(len(x)), m])
    rx = x - z @ np.linalg.lstsq(z, x, rcond=None)[0]
    ry = y - z @ np.linalg.lstsq(z, y, rcond=None)[0]
    return float(np.corrcoef(rx, ry)[0, 1])


def load_cohort(traits_zip: Path, env_zip: Path):
    endpoints = sorted({col for _, cols, _ in C.values() for col in cols})
    if len(endpoints) != 18:
        raise ValueError("Frozen complete-18 definition drift")
    with ZipFile(env_zip) as z:
        with z.open("process_environment/strict_spatial_chelsa_process.csv") as f:
            env = pd.read_csv(f, usecols=[
                "obs_id", "latitude", "longitude", "chelsa_bio12",
                "chelsa_rsds_mean"], dtype={"obs_id": str})
    chosen = set(env.obs_id)
    chunks = []
    with ZipFile(traits_zip) as z:
        with z.open("universe/continuous_trait_universe_observation_long.csv") as f:
            for x in pd.read_csv(f, chunksize=150000, low_memory=False, usecols=[
                    "obs_id", "taxon_name", "endpoint_id", "value"],
                    dtype={"obs_id": str, "taxon_name": str, "endpoint_id": str}):
                chunks.append(x[x.endpoint_id.isin(endpoints) &
                                x.obs_id.isin(chosen)])
        with z.open("extended_measurements/extended_continuous_head_level.csv") as f:
            photo = pd.read_csv(f, low_memory=False, usecols=[
                "obs_id", "taxon_name", "min_dimension", "sharpness",
                "mask_quality"], dtype={"obs_id": str, "taxon_name": str})
    t = pd.concat(chunks, ignore_index=True)
    med = t.groupby(["taxon_name", "endpoint_id"]).value.agg(
        ["median", "count"]).reset_index()
    eligible = med.loc[med["count"] >= 5]
    med = eligible.pivot(index="taxon_name", columns="endpoint_id",
                         values="median")
    means, scales = med.mean(), med.std(ddof=0)
    wide = t.pivot(index=["obs_id", "taxon_name"],
                   columns="endpoint_id", values="value").reset_index()
    wide = wide.dropna(subset=endpoints)
    counts = wide.groupby("taxon_name").size()
    wide = wide[wide.taxon_name.isin(counts[counts >= 5].index)].copy()
    wide = wide.sort_values(["taxon_name", "obs_id"]).reset_index(drop=True)
    if len(wide) != 1734 or wide.taxon_name.nunique() != 42:
        raise ValueError("Frozen 1734/42 common cohort did not reproduce")
    labels = wide.taxon_name.to_numpy()
    names = sorted(wide.taxon_name.unique())
    X = {}
    for name, (kind, columns, signs) in C.items():
        if kind == "scalar":
            if min(scales[columns]) <= 0:
                raise ValueError(f"Constant input scalar member {name}")
            score = ((wide[columns] - means[columns]) /
                     scales[columns]).to_numpy(float)
            X[name] = (score @ np.asarray(signs))[:, None] / len(columns)
        else:
            X[name] = wide[columns].to_numpy(float)
    covariates = wide[["obs_id", "taxon_name"]].merge(
        env, on="obs_id", validate="one_to_one").groupby(
            "taxon_name")[["latitude", "chelsa_bio12",
                            "chelsa_rsds_mean"]].median().loc[names]
    ph = photo.loc[photo.obs_id.isin(set(wide.obs_id))]
    if ph.obs_id.nunique() != len(wide):
        raise ValueError("Photo-quality sample did not match common cohort")
    extra = ph.groupby("taxon_name")[[
        "min_dimension", "sharpness", "mask_quality"]].median().loc[names]
    covariates = pd.concat([covariates, extra], axis=1)
    if not np.isfinite(covariates.to_numpy(float)).all():
        raise ValueError("Incomplete photo-quality sensitivity covariates")
    return X, labels, names, covariates


def reconstruct(xdict, labels, names, covariates, n_boot, n_perm, seed):
    n = len(names)
    rng = np.random.default_rng(seed)
    inds = rng.integers(0, n, size=(n_boot, n))
    orders = np.argsort(rng.random((n_perm, n)), axis=1)
    counts = pd.Series(labels).value_counts()
    weights = np.array([1 / counts[t] for t in labels], dtype=float)
    X = {k: pd.DataFrame(value).groupby(labels).median().loc[names].to_numpy()
         for k, value in xdict.items()}
    centered = {
        k: value - pd.DataFrame(value).groupby(labels).transform("mean").to_numpy()
        for k, value in xdict.items()
    }
    rows = []
    order = list(C)
    for i, left in enumerate(order):
        xa, xc = X[left], centered[left]
        pxx = np.stack([(xc[labels == t].T @ xc[labels == t]) /
                        sum(labels == t) for t in names])
        for right in order[i + 1:]:
            ya, yc = X[right], centered[right]
            pyy = np.stack([(yc[labels == t].T @ yc[labels == t]) /
                            sum(labels == t) for t in names])
            pxy = np.stack([(xc[labels == t].T @ yc[labels == t]) /
                            sum(labels == t) for t in names])
            w = rv(xc, yc, weights)
            a = rv(xa, ya)
            ax, by = xa[inds], ya[inds]
            dx, dy = ax - ax.mean(1, keepdims=True), by - by.mean(1, keepdims=True)
            a_boot = batch_rv(
                np.einsum("bti,btj->bij", dx, dx, optimize=True)/n,
                np.einsum("bti,btj->bij", dy, dy, optimize=True)/n,
                np.einsum("bti,btj->bij", dx, dy, optimize=True)/n)
            w_boot = batch_rv(pxx[inds].mean(1), pyy[inds].mean(1),
                              pxy[inds].mean(1))
            xx = xa - xa.mean(0)
            yy = ya - ya.mean(0)
            xy_perm = np.einsum(
                "ti,btj->bij", xx, yy[orders], optimize=True) / n
            den = np.sqrt(np.sum((xx.T@xx/n)**2) *
                          np.sum((yy.T@yy/n)**2))
            null = np.sum(xy_perm**2, axis=(-2, -1))/den
            pval = (1 + int(np.sum(null >= a - 1e-14))) / (n_perm + 1)
            jack = [rv(xa[np.arange(n) != k], ya[np.arange(n) != k])
                    for k in range(n)]
            rows.append({
                "construct_left": left, "construct_right": right,
                "within_rv": w, "among_rv": a, "delta_among_minus_within": a-w,
                "among_ci95_bootstrap": np.quantile(a_boot, [.025, .975]).tolist(),
                "within_ci95_bootstrap": np.quantile(w_boot, [.025, .975]).tolist(),
                "delta_ci95_bootstrap": np.quantile(
                    a_boot-w_boot, [.025, .975]).tolist(),
                "bootstrap_fraction_delta_positive": float(np.mean(a_boot > w_boot)),
                "permutation_p_independent_taxon_labels": pval,
                "jackknife_among_rv_min": float(min(jack)),
                "jackknife_among_rv_max": float(max(jack))
            })
    # Benjamini-Hochberg family across all 36 associations.
    ps = np.asarray([r["permutation_p_independent_taxon_labels"] for r in rows])
    idx = np.argsort(ps)
    q = np.minimum.accumulate((ps[idx] * len(ps) /
                               np.arange(1, len(ps)+1))[::-1])[::-1]
    out = np.empty(len(ps), dtype=float)
    out[idx] = np.minimum(q, 1)
    for r, qval in zip(rows, out):
        r["permutation_q_bh_36"] = float(qval)
    selected = []
    co = covariates[QC].to_numpy(float)
    for left, right in PRIMARY:
        item = next(r for r in rows if (r["construct_left"], r["construct_right"]) ==
                    (left, right))
        if X[left].shape[1] == X[right].shape[1] == 1:
            x, y = X[left].ravel(), X[right].ravel()
            c = partial_corr(x, y, co)
            signed = float(pearsonr(x, y).statistic)
            # Fresh, fixed taxon resampling for covariate-adjusted sensitivity.
            ri = np.random.default_rng(seed).integers(0, n, size=(4000, n))
            boot = [partial_corr(x[ii], y[ii], co[ii]) for ii in ri]
            boot = np.asarray(boot, dtype=float)
            item["signed_taxon_correlation"] = signed
            item["quality_climate_partial_r"] = c
            item["partial_r_bootstrap_ci95"] = np.quantile(
                boot[np.isfinite(boot)], [.025, .975]).tolist()
            item["partial_r_covariates"] = QC
        selected.append(item)
    return {
        "status": "DESCRIPTIVE_HEAD_COUPLING_PASS_NOT_ADAPTATION",
        "source_artifact_sha256": INPUTS,
        "n_taxa": n, "n_complete_observations": len(labels),
        "n_pairs": len(rows), "n_permutation_bh_36_lt05": sum(
            r["permutation_q_bh_36"] < .05 for r in rows),
        "n_among_rv_gt_within_rv": sum(
            r["among_rv"] > r["within_rv"] for r in rows),
        "bootstrap_draws": n_boot, "permutations_per_pair": n_perm,
        "selected_contrasts": selected, "all_pair_robustness": rows,
        "taxon_level_context_profiles": {
            "taxon_labels": names,
            "scalar_construct_medians": {
                construct: X[construct].ravel().tolist()
                for construct in ("presentation_angle", "head_elongation",
                                  "floral_chroma", "projection_prominence")
            },
            "matched_environment_and_image_quality": {
                field: covariates[field].to_numpy(float).tolist()
                for field in ("latitude", "longitude", "chelsa_bio12",
                              "chelsa_rsds_mean", "min_dimension",
                              "sharpness", "mask_quality")
            },
            "scope": "42 aggregate taxon medians, NOT species-specific biological interaction or native-range sampling."
        },
        "source_limits": [
            "Taxon-label null is not a lineage-history, phylogenetic correlation, or selection null.",
            "Taxon bootstrap assumes observed taxa as exchangeable; independent evolution is unverified.",
            "The angle is image-vertical and bract projection is silhouette roughness, NOT true gravitational inclination or actual spine length.",
            "Signed photo/geo residuals are post-hoc checks, not causal inference.",
            "No flower visitation effectiveness, oviposition, parasitoid early killing, or viable seed outcome in these inputs."
        ]
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--traits-zip", type=Path, required=True)
    p.add_argument("--environment-zip", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--bootstrap", type=int, default=3000)
    p.add_argument("--permutations", type=int, default=9999)
    p.add_argument("--seed", type=int, default=20261008)
    args = p.parse_args()
    if args.bootstrap < 100 or args.permutations < 99:
        p.error("Insufficient repetitions")
    verify_zip(args.traits_zip, "traits_zip")
    verify_zip(args.environment_zip, "environment_zip")
    X, labels, taxa, qc = load_cohort(
        args.traits_zip, args.environment_zip)
    result = reconstruct(X, labels, taxa, qc,
                         args.bootstrap, args.permutations, args.seed)
    assert result["n_pairs"] == 36
    assert result["n_complete_observations"] == 1734
    assert result["n_among_rv_gt_within_rv"] == 33
    assert result["n_permutation_bh_36_lt05"] == 7
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({
        "status": result["status"],
        "taxa": result["n_taxa"],
        "BH_significant_pairs": result["n_permutation_bh_36_lt05"],
        "selected": result["selected_contrasts"]
    }, indent=2))


if __name__ == "__main__":
    main()
