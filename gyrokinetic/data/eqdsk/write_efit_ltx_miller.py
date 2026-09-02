"""
Write an EFIT/GEQDSK equilibrium file for a Miller-geometry tokamak.

Geometry conventions follow rt_gk_ltx_aiwl_2x2v_p1.c (context-struct pattern):
all geometry functions take an explicit MillerParams argument instead of
capturing module-level globals.

Typical use
-----------
Run with defaults (reproduces the original LTX file):

    python write_efit_ltx_miller.py

Vary shape and resolution from the command line:

    python write_efit_ltx_miller.py --kappa 1.5 --delta 0.2 --nw 129 --nh 129 \
        --output shaped.geqdsk --plot

Supply an arbitrary q-profile from Python (programmatic use):

    from write_efit_ltx_miller import MillerParams, main
    params = MillerParams(qprofile_func=lambda r: 2.5 + 1.8 * r**2)
    main(params=params, output_file="custom_q.geqdsk")
"""

from __future__ import annotations

import argparse
import importlib
import math
from dataclasses import dataclass, field
from typing import Callable, Optional

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import scipy.optimize as sco
from scipy.interpolate import griddata

from eqdsk_io import write_geqdsk, write_wall_limiter


# ---------------------------------------------------------------------------
# Parameter dataclass
# ---------------------------------------------------------------------------

@dataclass
class MillerParams:
    """Physical parameters describing a Miller-geometry tokamak equilibrium.

    All input fields have defaults matching the LTX experiment.  Derived
    quantities are computed once in ``__post_init__`` and stored as read-only
    attributes.

    To use a fully custom safety-factor profile pass a callable to
    ``qprofile_func``; otherwise the analytic form parametrised by
    ``q_LCFS`` and ``s_LCFS`` is used.
    """

    # --- Primary inputs ---
    R_axis: float = 0.406051632
    """Major radius at the magnetic axis [m]."""

    B_axis: float = 0.240606108
    """Magnetic field magnitude at the axis [T]."""

    R_LCFSmid: float = 0.61183
    """Major radius of the LCFS at the outboard midplane [m]."""

    Rmid_min: Optional[float] = None
    """Minimum midplane major radius of the simulation box [m].
    Defaults to ``R_axis`` (axis at inner wall)."""

    Rmid_max: float = 0.7
    """Maximum midplane major radius of the simulation box [m]."""

    kappa: float = 1.3
    """Elongation (= 1 for circular cross-section)."""

    delta: float = 0.4
    """Triangularity (= 0 for no triangular shaping)."""

    q_LCFS: float = 4.3153848
    """Safety factor at the LCFS.  Used only when ``qprofile_func`` is None."""

    s_LCFS: float = 2.6899871
    """Magnetic shear at the LCFS.  Used only when ``qprofile_func`` is None."""

    qprofile_func: Optional[Callable[[float], float]] = field(
        default=None, repr=False
    )
    """Optional user-supplied q(r) callable.  Receives minor radius r [m] and
    returns the safety factor.  When provided, overrides the analytic form."""

    # --- Derived (computed in __post_init__) ---
    a_mid: float = field(init=False)
    R0: float = field(init=False)
    r0: float = field(init=False)
    B0: float = field(init=False)
    x_LCFS: float = field(init=False)
    q0: float = field(init=False)
    q_min: float = field(init=False)

    def __post_init__(self) -> None:
        if self.Rmid_min is None:
            self.Rmid_min = self.R_axis
        self.a_mid = self.R_LCFSmid - self.R_axis
        self.R0 = 0.5 * (self.Rmid_min + self.Rmid_max)
        self.r0 = self.R0 - self.R_axis
        self.B0 = self.B_axis * (self.R_axis / self.R0)
        self.x_LCFS = self.R_LCFSmid - self.Rmid_min
        self.q0 = self.qprofile(self.r0)
        self.q_min = self.qprofile(self.Rmid_min - self.R_axis)

    def qprofile(self, r: float) -> float:
        """Safety factor as a function of minor radius r [m]."""
        if self.qprofile_func is not None:
            return self.qprofile_func(r)
        return self.q_LCFS / (1.0 - self.s_LCFS * (r - self.a_mid) / self.a_mid)


# ---------------------------------------------------------------------------
# Geometry functions
# All follow the same convention as the C reference (rt_gk_ltx_aiwl_2x2v_p1.c):
# explicit params argument, no global captures.
# ---------------------------------------------------------------------------

def r_x(x: float, params: MillerParams) -> float:
    """Convert x coordinate (distance from Rmid_min) to minor radius r [m]."""
    return params.Rmid_min - params.R_axis + x


def R_f(r: float, theta: float, params: MillerParams) -> float:
    """Major radius R(r, theta) in the Miller parameterisation [m]."""
    return params.R_axis + r * np.cos(theta + np.arcsin(params.delta) * np.sin(theta))


def Z_f(r: float, theta: float, params: MillerParams) -> float:
    """Vertical position Z(r, theta) in the Miller parameterisation [m]."""
    return params.kappa * r * np.sin(theta)


def R_f_r(_r: float, theta: float, params: MillerParams) -> float:
    """Partial derivative dR/dr (independent of r in Miller geometry)."""
    return np.cos(theta + np.arcsin(params.delta) * np.sin(theta))


def R_f_theta(r: float, theta: float, params: MillerParams) -> float:
    """Partial derivative dR/dtheta."""
    asin_d = np.arcsin(params.delta)
    return -r * (asin_d * np.cos(theta) + 1.0) * np.sin(asin_d * np.sin(theta) + theta)


def Z_f_r(_r: float, theta: float, params: MillerParams) -> float:
    """Partial derivative dZ/dr (independent of r in Miller geometry)."""
    return params.kappa * np.sin(theta)


def Z_f_theta(r: float, theta: float, params: MillerParams) -> float:
    """Partial derivative dZ/dtheta."""
    return params.kappa * r * np.cos(theta)


def Jr_f(r: float, theta: float, params: MillerParams) -> float:
    """Jacobian determinant multiplied by R: R*(dR/dr * dZ/dtheta - dZ/dr * dR/dtheta)."""
    return R_f(r, theta, params) * (
        R_f_r(r, theta, params) * Z_f_theta(r, theta, params)
        - Z_f_r(r, theta, params) * R_f_theta(r, theta, params)
    )


def _integrand(t: float, r: float, params: MillerParams) -> float:
    """Integrand Jr/R² used in the flux calculations."""
    return Jr_f(r, t, params) / np.power(R_f(r, t, params), 2)


def dPsidr_f(r: float, _theta: float, params: MillerParams) -> float:
    """Radial flux gradient dPsi/dr at minor radius r [Wb/rad/m].

    The result is independent of theta because the integral runs over all
    poloidal angles [0, 2pi].  The theta argument is kept so that
    integrate.quad can pass it via args= when computing psi_fi.
    """
    integral, _ = integrate.quad(
        _integrand, 0.0, 2.0 * np.pi, args=(r, params), epsabs=1.0e-8
    )
    return params.B_axis * params.R_axis / (2.0 * np.pi * params.qprofile(r)) * integral


def psi_fi(r: float, theta: float, params: MillerParams) -> float:
    """Poloidal flux Psi(r) obtained by integrating dPsi/dr from 0 to r [Wb/rad]."""
    result, _ = integrate.quad(dPsidr_f, 0.0, r, args=(theta, params), epsabs=1.0e-8)
    return result


def Bphi_f(R: float, params: MillerParams) -> float:
    """Vacuum toroidal field B_phi = B_axis * R_axis / R [T]."""
    return params.B_axis * params.R_axis / R


def J_f(r: float, theta: float, params: MillerParams) -> float:
    """Full Jacobian Jr/dPsi_dr."""
    return Jr_f(r, theta, params) / dPsidr_f(r, theta, params)


# ---------------------------------------------------------------------------
# Field-line angle and shift functions
# (used by check_miller_shift.py and related analysis)
# ---------------------------------------------------------------------------

def alpha(r: float, theta: float, phi: float, params: MillerParams) -> float:
    """Field-line label angle alpha(r, theta, phi)."""
    twrap = (theta + math.pi) % (2 * math.pi) - math.pi

    def _wrap(t: float) -> float:
        return _integrand(t, r, params)

    if twrap > 0:
        res, _ = integrate.quad(_wrap, 0.0, twrap, epsabs=1e-10)
    else:
        res, _ = integrate.quad(_wrap, twrap, 0.0, epsabs=1e-10)
        res = -res

    return phi - (params.B_axis * params.R_axis * res) / dPsidr_f(r, theta, params)


def shift_lo_i(r: float, params: MillerParams) -> float:
    """Shift at the lower boundary."""
    return -params.r0 / params.q0 * alpha(r, -np.pi, 0.0, params)


def shift_lo_angle_i(r: float, params: MillerParams) -> float:
    """Field-line shift angle at the lower poloidal boundary."""
    return -alpha(r, -np.pi, 0.0, params)


def shift_up_angle_i(r: float, params: MillerParams) -> float:
    """Field-line shift angle at the upper poloidal boundary."""
    return -alpha(r, np.pi, 0.0, params)


def shift_lo_angle_i_basic(r: float, params: MillerParams) -> float:
    """Approximate shift angle (2*pi*q approximation)."""
    return -2 * np.pi * params.qprofile(r)


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

def main(
    params: Optional[MillerParams] = None,
    output_file: str = "ltx_miller.geqdsk",
    nw: int = 81,
    nh: int = 91,
    rmin: float = 0.1,
    rmax: float = 0.7,
    zmin: float = -0.4,
    zmax: float = 0.4,
    show_plots: bool = False,
    write_wall: bool = False,
) -> None:
    """Generate a Miller-geometry GEQDSK equilibrium file.

    Parameters
    ----------
    params      : Miller geometry parameters.  Defaults to LTX values.
    output_file : Path for the output GEQDSK file.
    nw          : Number of grid points in R.
    nh          : Number of grid points in Z.
    rmin, rmax  : R domain bounds [m].
    zmin, zmax  : Z domain bounds [m].
    show_plots  : Display diagnostic contour/scatter plots (blocking).
    write_wall  : Write companion wall and limiter ASCII files.
    """
    if params is None:
        params = MillerParams()

    # --- Grid setup ---
    rdim = rmax - rmin
    zdim = zmax - zmin
    rleft = rmin
    zmid = 0.5 * (zmin + zmax)

    xmin = 0.0
    xmax = params.Rmid_max - params.Rmid_min
    x = np.linspace(xmin, xmax, nw)
    r = params.Rmid_min - params.R_axis + x
    theta = np.linspace(-np.pi, np.pi, 65)
    r_grid, theta_grid = np.meshgrid(r, theta, indexing="ij")
    Rm = R_f(r_grid, theta_grid, params)
    Zm = Z_f(r_grid, theta_grid, params)

    # --- Compute psi(r) on 1-D radial grid ---
    psi = np.array([psi_fi(ri, 0.0, params) for ri in r])
    psi_plot = np.repeat(psi, len(theta)).reshape(len(r), len(theta))

    if show_plots:
        fig, ax = plt.subplots()
        cax = ax.contour(Rm, Zm, psi_plot, cmap="inferno")
        ax.set_title("Original Psi")
        plt.colorbar(cax)
        plt.show()

    # --- Interpolate to Cartesian (R, Z) grid ---
    Rgrid = np.linspace(rmin, rmax, nw)
    Zgrid = np.linspace(zmin, zmax, nh)
    grid_r, grid_z = np.meshgrid(Rgrid, Zgrid)
    flat_points = np.stack((Rm, Zm), axis=-1).reshape(-1, 2)
    flat_psi = psi_plot.flatten()
    psiRZ = griddata(flat_points, flat_psi, (grid_r, grid_z), method="cubic")
    nan_mask = np.isnan(psiRZ)
    psiRZ[nan_mask] = griddata(
        flat_points, flat_psi, (grid_r[nan_mask], grid_z[nan_mask]), method="nearest"
    )
    # psiRZ has shape (nh, nw); preserve transpose pair from original for PSIZR
    PSIZR = psiRZ.T.T  # (nh, nw) — no-op, kept explicit for clarity

    if show_plots:
        fig, ax = plt.subplots()
        cax = ax.contour(Rgrid, Zgrid, PSIZR, cmap="inferno")
        ax.set_title("Interpolated Psi")
        plt.colorbar(cax)
        plt.show()

    # --- Scalar equilibrium quantities ---
    SIMAG = psi_fi(0.0, 0.0, params)
    SIBRY = psi_fi(params.Rmid_min - params.R_axis + params.x_LCFS, 0.0, params)
    RMAXIS = params.R_axis
    ZMAXIS = 0.0
    RCENTR = params.R_axis
    BCENTR = params.B_axis
    CURRENT = 0.0

    # --- 1-D profiles on uniform psi grid ---
    NPSI = nw
    PSIGRID = np.linspace(SIMAG, SIBRY, NPSI)
    FPOL = np.full(NPSI, params.B_axis * params.R_axis)
    FFPRIM = np.zeros(NPSI)
    PPRIME = np.full(NPSI, -1e-6)
    PRES = integrate.cumulative_trapezoid(PPRIME, PSIGRID, initial=0)

    # --- Safety factor profile via flux inversion ---
    r_LCFS = r_x(params.x_LCFS, params)

    def _rlossfunc(r_val: float, psi0: float) -> float:
        return psi0 - psi_fi(r_val, 0.0, params)

    def _rfunc(psi: float) -> float:
        return sco.ridder(_rlossfunc, 0.0, r_LCFS, args=(psi,))

    QPSI = np.array([params.qprofile(_rfunc(p)) for p in PSIGRID])

    # --- LCFS boundary and limiter ---
    LIMITR = 10
    RBBBS = Rm[-1, :]
    ZBBBS = Zm[-1, :]
    RLIM = np.linspace(rleft, R_f(r_LCFS, np.pi, params), LIMITR)
    ZLIM = np.zeros(LIMITR)

    if show_plots:
        plt.figure()
        plt.scatter(RLIM, ZLIM, label="limiter")
        plt.scatter(RBBBS, ZBBBS, label="boundary")
        plt.legend()
        plt.show()

    # --- Write GEQDSK file ---
    write_geqdsk(
        filename=output_file,
        nw=nw,
        nh=nh,
        rdim=rdim,
        zdim=zdim,
        rcentr=RCENTR,
        rleft=rleft,
        zmid=zmid,
        rmaxis=RMAXIS,
        zmaxis=ZMAXIS,
        simag=SIMAG,
        sibry=SIBRY,
        bcentr=BCENTR,
        current=CURRENT,
        fpol=FPOL,
        pres=PRES,
        ffprim=FFPRIM,
        pprime=PPRIME,
        psizr=PSIZR,
        qpsi=QPSI,
        rbbbs=RBBBS,
        zbbbs=ZBBBS,
        rlim=RLIM,
        zlim=ZLIM,
    )
    print(f"Written: {output_file}")

    if write_wall:
        base = output_file.removesuffix(".geqdsk").removesuffix(".eqdsk")
        wall_path = base + "_wall"
        lim_path = base + "_limiter"
        write_wall_limiter(wall_path, lim_path, RBBBS, ZBBBS, RLIM, ZLIM)
        print(f"Written: {wall_path}, {lim_path}")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate a Miller-geometry EFIT/GEQDSK equilibrium file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    geo = p.add_argument_group("geometry")
    geo.add_argument("--r-axis", type=float, default=0.406051632,
                     metavar="M", help="Major radius at magnetic axis [m]")
    geo.add_argument("--b-axis", type=float, default=0.240606108,
                     metavar="T", help="Magnetic field at axis [T]")
    geo.add_argument("--r-lcfs-mid", type=float, default=0.61183,
                     metavar="M", help="Major radius of LCFS at outboard midplane [m]")
    geo.add_argument("--rmid-min", type=float, default=None,
                     metavar="M", help="Min midplane major radius [m] (default: R_axis)")
    geo.add_argument("--rmid-max", type=float, default=0.7,
                     metavar="M", help="Max midplane major radius [m]")
    geo.add_argument("--kappa", type=float, default=1.3,
                     help="Elongation (1 = circular)")
    geo.add_argument("--delta", type=float, default=0.4,
                     help="Triangularity (0 = no shaping)")
    geo.add_argument("--q-lcfs", type=float, default=4.3153848,
                     help="Safety factor at LCFS (analytic profile only)")
    geo.add_argument("--s-lcfs", type=float, default=2.6899871,
                     help="Magnetic shear at LCFS (analytic profile only)")
    geo.add_argument(
        "--qprofile-module", type=str, default=None, metavar="MODULE:FUNC",
        help=(
            "Import a custom q(r) callable, e.g. 'my_profiles:q_custom'. "
            "Overrides --q-lcfs and --s-lcfs for the profile shape."
        ),
    )

    grid = p.add_argument_group("grid")
    grid.add_argument("--nw", type=int, default=81, help="Radial grid points")
    grid.add_argument("--nh", type=int, default=91, help="Vertical grid points")
    grid.add_argument("--rmin", type=float, default=0.1, metavar="M",
                      help="Left R boundary [m]")
    grid.add_argument("--rmax", type=float, default=0.7, metavar="M",
                      help="Right R boundary [m]")
    grid.add_argument("--zmin", type=float, default=-0.4, metavar="M",
                      help="Bottom Z boundary [m]")
    grid.add_argument("--zmax", type=float, default=0.4, metavar="M",
                      help="Top Z boundary [m]")

    out = p.add_argument_group("output")
    out.add_argument("--output", type=str, default="ltx_miller.geqdsk",
                     help="Output GEQDSK file path")
    out.add_argument("--write-wall", action="store_true",
                     help="Write companion wall and limiter ASCII files")

    disp = p.add_argument_group("display")
    disp.add_argument("--plot", action="store_true",
                      help="Show diagnostic plots (blocking)")

    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()

    qprofile_func: Optional[Callable[[float], float]] = None
    if args.qprofile_module is not None:
        module_name, func_name = args.qprofile_module.rsplit(":", 1)
        mod = importlib.import_module(module_name)
        qprofile_func = getattr(mod, func_name)

    params = MillerParams(
        R_axis=args.r_axis,
        B_axis=args.b_axis,
        R_LCFSmid=args.r_lcfs_mid,
        Rmid_min=args.rmid_min,
        Rmid_max=args.rmid_max,
        kappa=args.kappa,
        delta=args.delta,
        q_LCFS=args.q_lcfs,
        s_LCFS=args.s_lcfs,
        qprofile_func=qprofile_func,
    )

    main(
        params=params,
        output_file=args.output,
        nw=args.nw,
        nh=args.nh,
        rmin=args.rmin,
        rmax=args.rmax,
        zmin=args.zmin,
        zmax=args.zmax,
        show_plots=args.plot,
        write_wall=args.write_wall,
    )
