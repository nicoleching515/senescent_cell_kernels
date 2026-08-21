"""
Synthetic tissue generator for the SASP spatial response kernel identifiability
study (Master Plan Section 22, Step 1.1-1.3).

Design goals, from the plan:
  * cells placed by a spatial point process with realistic density
    (median nearest-neighbour distance ~10-20 um; Section 3, Section 8 Test 1)
  * senders designated with CONTROLLED CLUSTERING, Poisson -> Thomas-like
    cluster process, clustering strength sweepable (Section 22 Step 1.1)
  * a KNOWN planted kernel r_i = mu_c + beta*exp(-d_i/lambda) + gamma'z_i + eps
    with d_i = distance to NEAREST sender (Section 6.1)
  * realistic nuisance (Section 22 Step 1.3):
      (a) spatially autocorrelated baseline, Gaussian random field, controllable
          correlation length
      (b) cell-type-dependent baselines
      (c) count noise (heteroscedastic in total counts)
      (d) a SENDER-DENSITY-CORRELATED CONFOUNDER -- "the key".  A latent niche
          field u(x) drives (i) where senders sit, (ii) local cell density,
          (iii) local cell-type composition and (iv) the baseline response.
          u is NEVER observed by the estimators.  This is what makes the
          neighbourhoods around senders genuinely different for reasons that
          have nothing to do with signalling.

Engineering rules (Section 18): scipy.spatial.cKDTree for ALL neighbour work.
No (n, n) distance matrix is ever built anywhere in this file.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, field
from typing import Dict, Any

import numpy as np
from scipy.spatial import cKDTree

# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------


@dataclass
class TissueConfig:
    # --- window and cell placement -----------------------------------------
    window_um: float = 2000.0
    # Matern-II hard-core thinning: proposal intensity and hard-core radius.
    # Calibrated (see calibrate_density.py) to give median NN distance ~13 um.
    prop_intensity: float = 0.020        # proposed points per um^2
    hardcore_um: float = 9.0             # minimum centre-to-centre distance
    density_mod: float = 0.55            # strength of u -> local density
    retain_base: float = 0.80            # mean retention after density modulation

    # --- cell types --------------------------------------------------------
    n_types: int = 4
    type_field_len_um: float = 180.0     # spatial scale of tissue "regions"
    type_field_amp: float = 1.1
    type_u_load: tuple = (0.9, -0.6, 0.3, -0.4)   # cell type composition <- u

    # --- senders -----------------------------------------------------------
    prevalence: float = 0.05             # fraction of cells that are senders
    clustering: float = 0.0              # kappa: 0 = Poisson (random senders)
    n_parents: int = 45                  # Thomas-process parents
    parent_sigma_um: float = 70.0        # Thomas cluster sd
    sender_type_pref: tuple = (0.7, 0.0, -0.5, 0.0)   # senescence is type-biased

    # --- the planted kernel ------------------------------------------------
    lambda_true_um: float = 30.0
    beta_true: float = 1.0
    kernel_family: str = "exponential"   # exponential|gaussian|powerlaw|step
    kernel_p: float = 2.0                # power-law exponent
    superposition: bool = False          # Section 6.3: sum over ALL senders
    superposition_trunc: float = 6.0     # truncate the sum at trunc*lambda

    # --- nuisance ----------------------------------------------------------
    autocorr_len_um: float = 30.0        # ell: correlation length of u and b
    conf_strength: float = 1.0           # scales BOTH u->sender and u->response
    psi_u_sender: float = 1.30           # u -> sender propensity (x conf_strength)
    eta_u_response: float = 0.70         # u -> response          (x conf_strength)
    base_grf_amp: float = 0.55           # amplitude of the *unconfounded* GRF b
    mu_types: tuple = (0.0, 0.40, -0.30, 0.15)
    gamma_density: float = 0.25          # response <- standardized local density
    gamma_counts: float = 0.30           # response <- standardized log counts
    sigma_eps: float = 0.80              # base noise sd at median counts
    counts_median: float = 300.0
    counts_log_sd: float = 0.55

    # --- numerics ----------------------------------------------------------
    grf_grid: int = 1024                 # FFT grid for the random fields
    knn_k: int = 20                      # k for kNN composition covariate
    density_radius_um: float = 50.0      # Section 8 Test 5 matching covariate

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# --------------------------------------------------------------------------
# Gaussian random fields (spectral synthesis on a torus, O(M^2 log M))
# --------------------------------------------------------------------------


def grf_grid(rng: np.random.Generator, window_um: float, corr_len_um: float,
             grid_n: int) -> np.ndarray:
    """Unit-variance GRF on a periodic grid with exponential covariance
    C(r) = exp(-r/ell).  Uses the Whittle-Matern (nu = 1/2, d = 2) spectrum
    S(k) proportional to (1/ell^2 + |k|^2)^(-3/2).

    Exponential covariance is chosen so that `ell` is an e-folding length and
    is therefore directly comparable to lambda_true of the exp kernel: the
    Figure 1 axis is "baseline autocorrelation length relative to lambda_true".
    """
    h = window_um / grid_n
    kx = 2.0 * np.pi * np.fft.fftfreq(grid_n, d=h)
    ky = 2.0 * np.pi * np.fft.rfftfreq(grid_n, d=h)
    k2 = kx[:, None] ** 2 + ky[None, :] ** 2
    amp = (k2 + 1.0 / corr_len_um ** 2) ** (-0.75)   # sqrt of the spectrum
    white = rng.standard_normal((grid_n, grid_n))
    fld = np.fft.irfft2(np.fft.rfft2(white) * amp, s=(grid_n, grid_n))
    fld -= fld.mean()
    sd = fld.std()
    if sd > 0:
        fld /= sd
    return fld


def sample_grid(fld: np.ndarray, window_um: float, coords: np.ndarray) -> np.ndarray:
    """Bilinear interpolation of a periodic grid field at scattered points."""
    grid_n = fld.shape[0]
    h = window_um / grid_n
    gx = coords[:, 0] / h
    gy = coords[:, 1] / h
    i0 = np.floor(gx).astype(np.int64)
    j0 = np.floor(gy).astype(np.int64)
    fx = gx - i0
    fy = gy - j0
    i0 %= grid_n
    j0 %= grid_n
    i1 = (i0 + 1) % grid_n
    j1 = (j0 + 1) % grid_n
    return (fld[i0, j0] * (1 - fx) * (1 - fy)
            + fld[i1, j0] * fx * (1 - fy)
            + fld[i0, j1] * (1 - fx) * fy
            + fld[i1, j1] * fx * fy)


# --------------------------------------------------------------------------
# point process
# --------------------------------------------------------------------------


def planted_kernel(d: np.ndarray, family: str, lam: float,
                   p: float = 2.0) -> np.ndarray:
    """The kernel actually planted in the simulation (Section 6.2 families)."""
    if family == "exponential":
        return np.exp(-d / lam)
    if family == "gaussian":
        return np.exp(-0.5 * (d / lam) ** 2)
    if family == "powerlaw":
        return (1.0 + d / lam) ** (-p)
    if family == "step":
        return (d < lam).astype(float)
    raise ValueError(family)


def superposition_signal(coords: np.ndarray, sender_coords: np.ndarray,
                         family: str, lam: float, p: float,
                         trunc_mult: float, window: float) -> np.ndarray:
    """sum_j K(||x_i - x_j||) over ALL senders (Section 6.3).

    Built from a TRUNCATED SPARSE distance matrix between the cell tree and the
    sender tree, never a dense (n_cells, n_senders) or (n, n) array.  For the
    step kernel the truncation radius is lambda itself, so the sum is exact;
    for the decaying families the tail beyond trunc_mult*lambda is negligible.
    """
    trunc = lam if family == "step" else trunc_mult * lam
    tc = cKDTree(coords, boxsize=window)
    ts = cKDTree(sender_coords, boxsize=window)
    D = tc.sparse_distance_matrix(ts, trunc, output_type="coo_matrix")
    vals = planted_kernel(D.data, family, lam, p)
    out = np.zeros(coords.shape[0])
    np.add.at(out, D.row, vals)
    return out


def hardcore_points(rng: np.random.Generator, window_um: float,
                    prop_intensity: float, hardcore_um: float) -> np.ndarray:
    """Matern type-II hard-core process: dense Poisson proposal, then delete any
    point that has a neighbour within `hardcore_um` carrying a smaller mark.
    Gives tissue-like packing (cells cannot overlap) rather than a raw Poisson
    process, which would produce unrealistically tiny nearest-neighbour
    distances.  Uses cKDTree.query_pairs -- never an (n, n) matrix.
    """
    n_prop = rng.poisson(prop_intensity * window_um ** 2)
    pts = rng.uniform(0.0, window_um, size=(n_prop, 2))
    marks = rng.random(n_prop)
    tree = cKDTree(pts, boxsize=window_um)          # periodic, avoids edge bias
    pairs = tree.query_pairs(hardcore_um, output_type='ndarray')
    keep = np.ones(n_prop, dtype=bool)
    if pairs.size:
        a, b = pairs[:, 0], pairs[:, 1]
        loser = np.where(marks[a] > marks[b], a, b)
        keep[loser] = False
    return pts[keep]


def _gumbel_topk(rng: np.random.Generator, log_w: np.ndarray, k: int) -> np.ndarray:
    """Weighted sampling without replacement (Plackett-Luce) via the Gumbel
    top-k trick.  Exact and O(n log n)."""
    g = rng.gumbel(size=log_w.shape[0])
    return np.argpartition(-(log_w + g), k - 1)[:k]


# --------------------------------------------------------------------------
# the generator
# --------------------------------------------------------------------------


def simulate_tissue(cfg: TissueConfig, seed) -> Dict[str, Any]:
    """Generate one synthetic tissue section.

    `seed` may be an int or a sequence of ints (passed to SeedSequence), so a
    sweep cell index and a replicate index can be combined reproducibly.
    """
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    W = cfg.window_um

    # ---- 1. cell positions -------------------------------------------------
    coords = hardcore_points(rng, W, cfg.prop_intensity, cfg.hardcore_um)

    # ---- 2. latent fields --------------------------------------------------
    ell = cfg.autocorr_len_um
    u_grid = grf_grid(rng, W, ell, cfg.grf_grid)          # the CONFOUNDER field
    b_grid = grf_grid(rng, W, ell, cfg.grf_grid)          # pure baseline GRF
    u = sample_grid(u_grid, W, coords)
    b = sample_grid(b_grid, W, coords)

    # ---- 3. density modulated by u  (local density becomes a real proxy for u)
    p_keep = cfg.retain_base * np.exp(cfg.density_mod * (u - u.mean()))
    p_keep = np.clip(p_keep, 0.02, 0.995)
    keep = rng.random(coords.shape[0]) < p_keep
    coords = coords[keep]
    u = u[keep]
    b = b[keep]
    n = coords.shape[0]

    # ---- 4. cell types, spatially patchy and correlated with u -------------
    logits = np.empty((n, cfg.n_types))
    for k in range(cfg.n_types):
        vk = sample_grid(grf_grid(rng, W, cfg.type_field_len_um, 256), W, coords)
        logits[:, k] = cfg.type_field_amp * vk + cfg.type_u_load[k] * u
    logits -= logits.max(axis=1, keepdims=True)
    p = np.exp(logits)
    p /= p.sum(axis=1, keepdims=True)
    cum = np.cumsum(p, axis=1)
    cell_type = (rng.random((n, 1)) > cum).sum(axis=1).astype(np.int8)

    # ---- 5. senders with controlled clustering ----------------------------
    #   kappa = 0                -> senders are a random (Poisson) subset
    #   kappa > 0                -> Thomas-like aggregation around parents
    if cfg.clustering > 0 and cfg.n_parents > 0:
        parents = rng.uniform(0.0, W, size=(cfg.n_parents, 2))
        # Gaussian cluster kernel summed over parents, on the torus.
        # (n, n_parents) with n_parents ~ 45 -- small and vectorized, this is
        # NOT an (n, n) matrix.
        delta = np.abs(coords[:, None, :] - parents[None, :, :])
        delta = np.minimum(delta, W - delta)
        dd2 = (delta ** 2).sum(axis=2)
        g = np.exp(-0.5 * dd2 / cfg.parent_sigma_um ** 2).sum(axis=1)
        g = np.log1p(g)
        g = (g - g.mean()) / (g.std() + 1e-12)
    else:
        g = np.zeros(n)

    psi = cfg.psi_u_sender * cfg.conf_strength
    log_w = (cfg.clustering * g
             + psi * u
             + np.asarray(cfg.sender_type_pref)[cell_type])
    n_send = max(5, int(round(cfg.prevalence * n)))
    n_send = min(n_send, n - 10)
    sender_idx = _gumbel_topk(rng, log_w, n_send)
    sender_mask = np.zeros(n, dtype=bool)
    sender_mask[sender_idx] = True

    # ---- 6. geometry: KD-tree only (Section 18.1) --------------------------
    tree_all = cKDTree(coords)
    tree_send = cKDTree(coords[sender_mask])
    d_sender, _ = tree_send.query(coords, k=1, workers=1)

    dens50 = np.array(tree_all.query_ball_point(
        coords, r=cfg.density_radius_um, return_length=True)) - 1.0

    kk = min(cfg.knn_k + 1, n)
    _, knn_idx = tree_all.query(coords, k=kk, workers=1)
    knn_types = cell_type[knn_idx[:, 1:]]
    knn_comp = np.stack(
        [(knn_types == k).mean(axis=1) for k in range(cfg.n_types)], axis=1)

    d_edge = np.minimum.reduce([coords[:, 0], coords[:, 1],
                                W - coords[:, 0], W - coords[:, 1]])

    # ---- 7. counts and the response ---------------------------------------
    counts = cfg.counts_median * np.exp(
        cfg.counts_log_sd * rng.standard_normal(n))
    log_counts_z = (np.log(counts) - np.log(cfg.counts_median)) / cfg.counts_log_sd
    dens_z = (dens50 - dens50.mean()) / (dens50.std() + 1e-12)

    eta = cfg.eta_u_response * cfg.conf_strength
    if cfg.superposition:
        raw = superposition_signal(coords, coords[sender_mask],
                                   cfg.kernel_family, cfg.lambda_true_um,
                                   cfg.kernel_p, cfg.superposition_trunc, W)
    else:
        raw = planted_kernel(d_sender, cfg.kernel_family,
                             cfg.lambda_true_um, cfg.kernel_p)
    signal = cfg.beta_true * raw
    noise_sd = cfg.sigma_eps * np.sqrt(cfg.counts_median / counts)

    r = (np.asarray(cfg.mu_types)[cell_type]
         + signal
         + eta * u
         + cfg.base_grf_amp * b
         + cfg.gamma_density * dens_z
         + cfg.gamma_counts * log_counts_z
         + noise_sd * rng.standard_normal(n))

    # ---- 8. diagnostics ----------------------------------------------------
    nn_d, _ = tree_all.query(coords, k=2, workers=1)
    median_nn = float(np.median(nn_d[:, 1]))

    return dict(
        coords=coords, cell_type=cell_type, sender_mask=sender_mask,
        d_sender=d_sender, r=r,
        dens50=dens50, knn_comp=knn_comp, log_counts_z=log_counts_z,
        d_edge=d_edge, counts=counts,
        u_latent=u, b_latent=b, signal_true=signal, signal_raw=raw,
        n_cells=n, n_senders=int(sender_mask.sum()),
        median_nn_um=median_nn,
        prevalence_real=float(sender_mask.mean()),
        window_um=W, cfg=cfg,
    )


# --------------------------------------------------------------------------
# realized sender clustering statistic (Section 8 Test 4)
# --------------------------------------------------------------------------


def ripley_ratio(coords: np.ndarray, sender_mask: np.ndarray,
                 radius: float = 50.0) -> float:
    """Ratio of the observed number of sender-sender pairs within `radius` to
    the number expected under a random relabelling of cells.  1.0 = Poisson,
    >1 = aggregated.  This is the same statistic Section 8 Test 4(b) asks for
    on real data, so real tissue can be located on the Figure 1 axis.
    """
    tree_all = cKDTree(coords)
    tree_s = cKDTree(coords[sender_mask])
    obs = tree_s.count_neighbors(tree_s, radius) - sender_mask.sum()
    tot = tree_all.count_neighbors(tree_all, radius) - coords.shape[0]
    ns, nt = int(sender_mask.sum()), coords.shape[0]
    exp = tot * ns * (ns - 1) / (nt * (nt - 1))
    return float(obs / exp) if exp > 0 else np.nan
