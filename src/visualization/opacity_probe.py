"""Load and plot SKIRT OpacityProbe optical-depth cubes."""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil
from pathlib import Path
from typing import Sequence

import numpy as np
from astropy.io import fits


AU_PER_PC = 206264.80624709636
DEFAULT_SLICE_WAVELENGTHS_MICRON: tuple[float, ...] = (0.1, 0.55, 2.0, 10.0, 100.0, 500.0)


@dataclass(frozen=True)
class OpacityCube:
    """Container for a projected SKIRT OpacityProbe FITS cube.

    The data array is ordered as ``(n_lambda, ny, nx)`` after Astropy reads the
    FITS file. Each wavelength plane stores optical depth, ``tau_lambda``.
    """

    data: np.ndarray
    wavelength_micron: np.ndarray
    x: np.ndarray
    y: np.ndarray
    data_unit: str
    spatial_unit: str
    fits_path: Path
    header: fits.Header


def load_opacity_cube(
    fits_path: str | Path,
    spatial_unit: str = "pc",
    memmap: bool = True,
) -> OpacityCube:
    """Load a SKIRT ``OpacityProbe`` FITS cube.

    Parameters
    ----------
    fits_path:
        Path to a SKIRT ``*_tau.fits`` file.
    spatial_unit:
        Unit for the returned x and y coordinate arrays. The current helper
        supports ``"pc"`` and the original FITS unit.
    memmap:
        Keep large FITS cubes memory-mapped when possible.
    """

    fits_path = Path(fits_path)

    with fits.open(fits_path, memmap=memmap) as hdul:
        data = hdul[0].data
        if data is None:
            raise ValueError(f"No primary image data found in {fits_path}")

        data = np.asarray(data)
        if not memmap:
            data = np.array(data, copy=True)
        if data.ndim == 2:
            data = data[np.newaxis, :, :]
        if data.ndim != 3:
            raise ValueError(f"Expected a 2D image or 3D cube in {fits_path}; got shape {data.shape}")

        header = hdul[0].header.copy()
        wavelength = _read_wavelength_axis(hdul)

    if len(wavelength) != data.shape[0]:
        raise ValueError(
            f"Wavelength axis has {len(wavelength)} values, but FITS cube has {data.shape[0]} wavelength planes"
        )

    x = _axis_coordinates(header, axis=1, n=data.shape[2], target_unit=spatial_unit)
    y = _axis_coordinates(header, axis=2, n=data.shape[1], target_unit=spatial_unit)

    return OpacityCube(
        data=data,
        wavelength_micron=np.asarray(wavelength, dtype=float),
        x=x,
        y=y,
        data_unit=str(header.get("BUNIT", "")).strip(),
        spatial_unit=spatial_unit,
        fits_path=fits_path,
        header=header,
    )


def nearest_wavelength_index(cube: OpacityCube, wavelength_micron: float) -> int:
    """Return the index of the wavelength plane closest to ``wavelength_micron``."""

    if cube.wavelength_micron.size == 0:
        raise ValueError("The cube has no wavelength axis")
    return int(np.nanargmin(np.abs(cube.wavelength_micron - wavelength_micron)))


def plot_tau_slice(
    cube: OpacityCube,
    wavelength_micron: float,
    save_path: str | Path | None = None,
    cmap: str = "magma",
    show_tau_one: bool = True,
):
    """Plot one optical-depth map at the nearest available wavelength."""

    index = nearest_wavelength_index(cube, wavelength_micron)
    frame = np.asarray(cube.data[index], dtype=float)

    fig, ax = _plt().subplots(figsize=(5.2, 4.4), constrained_layout=True)
    image = _imshow(ax, frame, cube, cmap=cmap, norm=_log_norm(frame))
    if show_tau_one and np.nanmin(frame) <= 1.0 <= np.nanmax(frame):
        ax.contour(cube.x, cube.y, frame, levels=[1.0], colors="white", linewidths=0.8)
    ax.set_title(f"tau at {cube.wavelength_micron[index]:.3g} micron")
    _label_spatial_axes(ax, cube)
    fig.colorbar(image, ax=ax, label="Optical depth tau")
    _save(fig, save_path)
    return fig, ax


def plot_tau_slices(
    cube: OpacityCube,
    wavelengths_micron: Sequence[float] | None = None,
    save_path: str | Path | None = None,
    cmap: str = "magma",
):
    """Plot several optical-depth maps in one figure."""

    if wavelengths_micron is None:
        wavelengths_micron = _default_slice_wavelengths(cube)

    indices = [nearest_wavelength_index(cube, wavelength) for wavelength in wavelengths_micron]
    frames = [np.asarray(cube.data[index], dtype=float) for index in indices]

    fig, axes = _make_axes(len(frames), figsize_per_panel=(4.2, 3.8))
    norm = _shared_log_norm(frames)
    for ax, frame, index in zip(_iter_axes(axes), frames, indices):
        image = _imshow(ax, frame, cube, cmap=cmap, norm=norm)
        if np.nanmin(frame) <= 1.0 <= np.nanmax(frame):
            ax.contour(cube.x, cube.y, frame, levels=[1.0], colors="white", linewidths=0.7)
        ax.set_title(f"{cube.wavelength_micron[index]:.3g} micron")
        _label_spatial_axes(ax, cube)
        fig.colorbar(image, ax=ax, label="Optical depth tau")
    _hide_unused_axes(axes, len(frames))
    _save(fig, save_path)
    return fig, axes


def plot_transmission_slice(
    cube: OpacityCube,
    wavelength_micron: float,
    save_path: str | Path | None = None,
    cmap: str = "viridis",
):
    """Plot transmission, ``exp(-tau)``, at one wavelength."""

    index = nearest_wavelength_index(cube, wavelength_micron)
    transmission = np.exp(-np.asarray(cube.data[index], dtype=float))

    fig, ax = _plt().subplots(figsize=(5.2, 4.4), constrained_layout=True)
    image = _imshow(ax, transmission, cube, cmap=cmap)
    ax.set_title(f"transmission at {cube.wavelength_micron[index]:.3g} micron")
    _label_spatial_axes(ax, cube)
    fig.colorbar(image, ax=ax, label="Transmission exp(-tau)")
    _save(fig, save_path)
    return fig, ax


def plot_opaque_fraction(
    cube: OpacityCube,
    thresholds: Sequence[float] = (1.0, 10.0, 100.0),
    save_path: str | Path | None = None,
):
    """Plot the fraction of pixels with optical depth above each threshold."""

    fig, ax = _plt().subplots(figsize=(6.0, 4.2), constrained_layout=True)
    for threshold in thresholds:
        fraction = np.nanmean(cube.data > threshold, axis=(1, 2))
        ax.plot(cube.wavelength_micron, fraction, label=f"tau > {threshold:g}")

    ax.set_xscale("log")
    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel("Fraction of pixels")
    ax.set_ylim(0.0, 1.0)
    ax.legend()
    _save(fig, save_path)
    return fig, ax


def plot_tau_percentiles(
    cube: OpacityCube,
    percentiles: Sequence[float] = (50.0, 90.0, 99.0),
    save_path: str | Path | None = None,
):
    """Plot optical-depth percentiles as a function of wavelength."""

    fig, ax = _plt().subplots(figsize=(6.0, 4.2), constrained_layout=True)
    for percentile in percentiles:
        tau = np.nanpercentile(cube.data, percentile, axis=(1, 2))
        ax.plot(cube.wavelength_micron, tau, label=f"{percentile:g}th percentile")

    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1e-3)
    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel("Optical depth tau")
    ax.legend()
    _save(fig, save_path)
    return fig, ax


def _read_wavelength_axis(hdul: fits.HDUList) -> np.ndarray:
    for hdu in hdul[1:]:
        if hdu.data is None:
            continue
        names = getattr(hdu.data, "names", None) or []
        if "GRID_POINTS" in names:
            return np.asarray(hdu.data["GRID_POINTS"], dtype=float)
        if names:
            return np.asarray(hdu.data[names[0]], dtype=float)
    raise ValueError("No wavelength table extension with GRID_POINTS was found")


def _axis_coordinates(header: fits.Header, axis: int, n: int, target_unit: str) -> np.ndarray:
    crpix = float(header.get(f"CRPIX{axis}", (n + 1) / 2.0))
    crval = float(header.get(f"CRVAL{axis}", 0.0))
    cdelt = float(header.get(f"CDELT{axis}", 1.0))
    source_unit = str(header.get(f"CUNIT{axis}", target_unit)).strip()
    coords = (np.arange(n, dtype=float) + 1.0 - crpix) * cdelt + crval

    if source_unit.lower() == target_unit.lower():
        return coords
    if source_unit.lower() == "au" and target_unit.lower() == "pc":
        return _convert_au_to_pc(coords)
    raise ValueError(f"Cannot convert spatial axis from {source_unit!r} to {target_unit!r}")


def _convert_au_to_pc(values: np.ndarray) -> np.ndarray:
    return values / AU_PER_PC


def _default_slice_wavelengths(cube: OpacityCube) -> list[float]:
    in_range = [
        wavelength
        for wavelength in DEFAULT_SLICE_WAVELENGTHS_MICRON
        if cube.wavelength_micron.min() <= wavelength <= cube.wavelength_micron.max()
    ]
    if in_range:
        return in_range
    return [float(np.nanmedian(cube.wavelength_micron))]


def _make_axes(n: int, figsize_per_panel: tuple[float, float]):
    ncols = min(3, max(1, n))
    nrows = ceil(n / ncols)
    fig, axes = _plt().subplots(
        nrows,
        ncols,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows),
        squeeze=False,
        constrained_layout=True,
    )
    return fig, axes


def _iter_axes(axes):
    return np.asarray(axes, dtype=object).ravel()


def _hide_unused_axes(axes, used: int) -> None:
    for ax in _iter_axes(axes)[used:]:
        ax.set_visible(False)


def _imshow(ax, frame: np.ndarray, cube: OpacityCube, cmap: str, norm=None):
    extent = [cube.x.min(), cube.x.max(), cube.y.min(), cube.y.max()]
    return ax.imshow(frame, origin="lower", extent=extent, aspect="equal", cmap=cmap, norm=norm)


def _log_norm(frame: np.ndarray):
    positive = np.asarray(frame)[np.asarray(frame) > 0]
    if positive.size == 0:
        return None
    return _colors().LogNorm(vmin=float(np.nanmin(positive)), vmax=float(np.nanmax(positive)))


def _shared_log_norm(frames: Sequence[np.ndarray]):
    positive_values = []
    for frame in frames:
        values = np.asarray(frame)
        positive_values.append(values[values > 0])
    positive_values = [values for values in positive_values if values.size]
    if not positive_values:
        return None
    positive = np.concatenate(positive_values)
    return _colors().LogNorm(vmin=float(np.nanmin(positive)), vmax=float(np.nanmax(positive)))


def _label_spatial_axes(ax, cube: OpacityCube) -> None:
    ax.set_xlabel(f"x [{cube.spatial_unit}]")
    ax.set_ylabel(f"y [{cube.spatial_unit}]")


def _save(fig, save_path: str | Path | None) -> None:
    if save_path is None:
        return
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=200, bbox_inches="tight")


def _plt():
    import matplotlib.pyplot as plt

    return plt


def _colors():
    from matplotlib import colors

    return colors
