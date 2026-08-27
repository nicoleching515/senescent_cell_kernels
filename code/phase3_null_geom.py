#!/usr/bin/env python3
"""Phase 7 / correction C1 — in-tissue coordinate nulls for Phase 3.

WHY THIS FILE EXISTS
--------------------
`run_phase3_nulls.py` shifts and rotates the sender set over the whole-section
BOUNDING BOX:

    lo = coords.min(0); hi = coords.max(0); span = hi - lo
    lo + (pts - lo + rng.uniform(0, 1, 2) * span) % span

A liver section is not a rectangle.  Only 62-86 % of each section's bounding box
is occupied by tissue, so a whole-section wrap throws a large fraction of the
shifted senders into empty space.  `reports/CS_PHASE4.md` §2.4 says this in as
many words ("A torus shift on a whole section throws ~20 % of shifted cells into
the void outside the tissue") and Phase 4 fixed it by running everything on
solid-tissue tiles.  Phase 3 was never re-run.  This module supplies the
corrected geometry; `run_phase3_nulls.py --stage perm_c1` runs it.

THE VARIANTS (Phase 7 §2)
-------------------------
translation family
  N3_orig    whole-section bounding-box wrap -- the ORIGINAL, kept as reference
  N3_tile    wrap inside Phase-4-style solid-tissue tiles (`phase4_tiles.tiles_for`)
  N3_occ     whole-section wrap, offset REJECTED unless >= 95 % of shifted
             senders land in an occupied 25 um grid cell   (the §2 spec)
  N3_occ15   ditto at >= 85 %                              (supplementary: the
             5 % criterion turns out to admit only near-identity offsets)
  N3_swap    every sender relocated to a randomly chosen REAL cell position
  N3_snap    whole-section wrap, then every shifted sender snapped to the
             nearest real cell position                    (supplementary)
rotation family -- identical constructions about the centroid
  N4_orig, N4_tile, N4_occ, N4_occ15, N4_swap (= rotate then snap), N4_snap

`N3_swap` relocates senders to uniformly chosen real cell positions, so unlike
N3 it does NOT preserve sender clustering; it is in-tissue by construction but
it is a strictly stronger null.  The literal swap is orientation-free, so there
is no distinct rotation analogue: `N4_swap` is defined as rotate-then-snap,
which keeps the rotated configuration and is in tissue by construction.

Every construction below returns a (k, 2) array of shifted SENDER coordinates.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

from phase4_tiles import tiles_for, TILE_UM      # Phase 4's own tile chooser

GRID_UM = 25.0          # the occupancy grid of Phase 7 §2 / master plan Tier D
OCC_TOL = 0.05          # §2: reject if > 5 % of senders land outside
OCC_TOL_RELAXED = 0.15  # supplementary
TILE_SOLID = 0.98       # a tile is "solid tissue" if >= 98 % of it is tissue
N_ANGLE = 720           # rotation-angle grid for the occupancy screen

TRANSLATION = ("N3_orig", "N3_tile", "N3_occ", "N3_occ15", "N3_swap", "N3_snap")
ROTATION = ("N4_orig", "N4_tile", "N4_occ", "N4_occ15", "N4_swap")
ALL_NULLS = TRANSLATION + ROTATION
TILE_NULLS = ("N3_tile", "N4_tile")
FULL_NULLS = tuple(n for n in ALL_NULLS if n not in TILE_NULLS)


class _Shim:
    """Minimal stand-in for phase4_data.Sec4 so `tiles_for` can be reused
    verbatim on a Phase 3 section."""

    def __init__(self, coords):
        self.coords = coords
        self.lo = coords.min(0)
        self.hi = coords.max(0)


class Geom:
    """All section geometry the corrected nulls need, built once per section."""

    def __init__(self, coords, sender, eligible, tile_side=TILE_UM,
                 grid_um=GRID_UM):
        self.coords = np.asarray(coords, float)
        self.sender = np.asarray(sender, bool)
        self.eligible = np.asarray(eligible, bool)
        self.n = self.coords.shape[0]
        self.k = int(self.sender.sum())
        self.grid_um = float(grid_um)
        self.lo = self.coords.min(0)
        self.hi = self.coords.max(0)
        self.span = self.hi - self.lo
        self.cen = self.coords.mean(0)
        self.tree = cKDTree(self.coords)
        self.elig_idx = np.flatnonzero(self.eligible)

        # ---- 25 um occupancy grid and the filled tissue mask --------------
        g = np.floor((self.coords - self.lo) / self.grid_um).astype(np.int64)
        self.nx, self.ny = int(g[:, 0].max()) + 1, int(g[:, 1].max()) + 1
        self.occ = np.zeros((self.nx, self.ny), bool)
        self.occ[g[:, 0], g[:, 1]] = True
        # tissue = occupancy closed and hole-filled: a 25 um cell with no cell
        # centroid in it but tissue all around it (a sinusoid, a vessel lumen)
        # is NOT the void outside the section.
        self.tissue = ndi.binary_fill_holes(
            ndi.binary_closing(self.occ, np.ones((5, 5), bool)))
        self.occ_frac = float(self.occ.mean())
        self.tissue_frac = float(self.tissue.mean())
        self._g = g
        self._frac = (self.coords - self.lo) / self.grid_um - g   # sub-cell part
        self.sg = g[self.sender]
        self.sfrac = self._frac[self.sender]

        # ---- accepted translation offsets, exactly, by FFT ---------------
        # corr[u] = fraction of senders that land in an OCCUPIED 25 um cell when
        # the sender grid is wrapped by the integer offset u.  This is a
        # circular cross-correlation, so one FFT gives the acceptance map over
        # every one of nx*ny candidate offsets -- no rejection loop needed.
        S = np.zeros((self.nx, self.ny), float)
        np.add.at(S, (self.sg[:, 0], self.sg[:, 1]), 1.0)
        corr = np.fft.irfft2(np.conj(np.fft.rfft2(S)) * np.fft.rfft2(
            self.occ.astype(float)), s=(self.nx, self.ny)) / max(self.k, 1)
        self.offset_in_occ = np.clip(corr, 0.0, 1.0)
        du = np.minimum(np.arange(self.nx), self.nx - np.arange(self.nx))
        dv = np.minimum(np.arange(self.ny), self.ny - np.arange(self.ny))
        self.offset_dist = np.hypot(du[:, None] * self.grid_um,
                                    dv[None, :] * self.grid_um)

        # ---- accepted rotation angles, by direct scan --------------------
        th = np.linspace(0.0, 2 * np.pi, N_ANGLE, endpoint=False)
        self.angles = th
        fr = np.empty(N_ANGLE)
        for i, t in enumerate(th):
            fr[i] = self._in_occ(self._rot_bbox_pts(t))
        self.angle_in_occ = fr

        # ---- Phase 4 solid-tissue tiles ----------------------------------
        self.tile_side = float(tile_side)
        self.tiles = self._solid_tiles()
        self.tile_of = np.full(self.n, -1, np.int32)
        for ti, (x0, y0) in enumerate(self.tiles):
            m = ((self.coords[:, 0] >= x0) & (self.coords[:, 0] < x0 + tile_side)
                 & (self.coords[:, 1] >= y0) & (self.coords[:, 1] < y0 + tile_side))
            self.tile_of[m] = ti
        self.in_tile = self.tile_of >= 0
        self.tile_cov_cells = float(self.in_tile.mean())
        self.tile_cov_senders = (float(self.in_tile[self.sender].mean())
                                 if self.k else np.nan)

    # ------------------------------------------------------------------ #
    def _solid_tiles(self):
        """Phase-4 tiles (`tiles_for`, side 1200 um, densest non-overlapping)
        kept only when the tile is >= TILE_SOLID tissue."""
        cand = tiles_for(_Shim(self.coords), n_tiles=10 ** 6,
                         side=self.tile_side)
        k = int(round(self.tile_side / self.grid_um))
        keep = []
        for _n, x0, y0 in cand:
            i0 = int(round((x0 - self.lo[0]) / self.grid_um))
            j0 = int(round((y0 - self.lo[1]) / self.grid_um))
            sub = self.tissue[i0:i0 + k, j0:j0 + k]
            if sub.shape == (k, k) and sub.mean() >= TILE_SOLID:
                keep.append((float(x0), float(y0)))
        return keep

    def _in_occ(self, pts):
        """Fraction of `pts` sitting in an occupied 25 um grid cell."""
        ij = np.floor((pts - self.lo) / self.grid_um).astype(np.int64)
        ij[:, 0] = np.clip(ij[:, 0], 0, self.nx - 1)
        ij[:, 1] = np.clip(ij[:, 1], 0, self.ny - 1)
        return float(self.occ[ij[:, 0], ij[:, 1]].mean())

    # ---- the ORIGINAL operations, copied verbatim from run_phase3_nulls --
    def _shift_bbox_pts(self, v01):
        return self.lo + (self.coords[self.sender] - self.lo
                          + v01 * self.span) % self.span

    def _rot_bbox_pts(self, th):
        R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        pts = self.coords[self.sender]
        return self.lo + ((pts - self.cen) @ R.T + self.cen - self.lo) % self.span

    # ------------------------------------------------------------------ #
    def accepted_offsets(self, tol):
        """Integer grid offsets whose in-occupancy fraction is >= 1 - tol."""
        A = np.argwhere(self.offset_in_occ >= 1.0 - tol)
        if A.size == 0:                     # nothing clears the bar
            A = np.argwhere(self.offset_in_occ >= self.offset_in_occ.max() - 1e-12)
        return A

    def accepted_angles(self, tol):
        ix = np.flatnonzero(self.angle_in_occ >= 1.0 - tol)
        if ix.size == 0:
            ix = np.array([int(np.argmax(self.angle_in_occ))])
        return ix

    # ---- null constructions ------------------------------------------- #
    def N3_orig(self, rng):
        return self._shift_bbox_pts(rng.uniform(0, 1, 2))

    def N4_orig(self, rng):
        return self._rot_bbox_pts(rng.uniform(0, 2 * np.pi))

    def _shift_lattice(self, u):
        """Wrap the sender grid by integer offset u, keeping the sub-cell part.
        Consistent with the FFT acceptance map by construction."""
        gg = np.empty_like(self.sg)
        gg[:, 0] = (self.sg[:, 0] + u[0]) % self.nx
        gg[:, 1] = (self.sg[:, 1] + u[1]) % self.ny
        return self.lo + (gg + self.sfrac) * self.grid_um

    def N3_occ(self, rng, tol=OCC_TOL):
        A = self._acc_off(tol)
        return self._shift_lattice(A[rng.integers(len(A))])

    def N3_occ15(self, rng):
        return self.N3_occ(rng, OCC_TOL_RELAXED)

    def N4_occ(self, rng, tol=OCC_TOL):
        ix = self._acc_ang(tol)
        return self._rot_bbox_pts(self.angles[ix[rng.integers(len(ix))]])

    def N4_occ15(self, rng):
        return self.N4_occ(rng, OCC_TOL_RELAXED)

    def _acc_off(self, tol):
        if not hasattr(self, "_ao"):
            self._ao = {}
        if tol not in self._ao:
            self._ao[tol] = self.accepted_offsets(tol)
        return self._ao[tol]

    def _acc_ang(self, tol):
        if not hasattr(self, "_aa"):
            self._aa = {}
        if tol not in self._aa:
            self._aa[tol] = self.accepted_angles(tol)
        return self._aa[tol]

    def N3_swap(self, rng):
        """Every sender relocated to a randomly chosen REAL cell position
        (drawn without replacement from the sender-eligible cells).  In tissue
        by construction; sender count preserved."""
        pick = rng.choice(self.elig_idx, size=self.k, replace=False)
        return self.coords[pick]

    def _snap(self, pts):
        _, ix = self.tree.query(pts, k=1, workers=1)
        return self.coords[ix]

    def N3_snap(self, rng):
        return self._snap(self.N3_orig(rng))

    def N4_swap(self, rng):
        """Rotation analogue of the swap: rotate, then snap every sender to the
        nearest real cell position."""
        return self._snap(self.N4_orig(rng))

    # ---- tile family ---------------------------------------------------- #
    def tile_sender_idx(self):
        return np.flatnonzero(self.sender & self.in_tile)

    def N3_tile(self, rng):
        """One shared random offset, wrapped INSIDE each solid tile."""
        v = rng.uniform(0, 1, 2) * self.tile_side
        ii = self.tile_sender_idx()
        org = np.array([self.tiles[t] for t in self.tile_of[ii]], float)
        return org + (self.coords[ii] - org + v) % self.tile_side

    def N4_tile(self, rng):
        """One shared random angle, about each solid tile's own centre,
        wrapped inside that tile."""
        th = rng.uniform(0, 2 * np.pi)
        R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        ii = self.tile_sender_idx()
        org = np.array([self.tiles[t] for t in self.tile_of[ii]], float)
        c = org + 0.5 * self.tile_side
        return org + (((self.coords[ii] - c) @ R.T + c - org) % self.tile_side)

    # ------------------------------------------------------------------ #
    def draw(self, null, rng):
        return getattr(self, null)(rng)
