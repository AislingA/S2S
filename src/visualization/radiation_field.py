"""Load and plot SKIRT RadiationFieldProbe output cubes."""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from astropy.io import fits


LIGHT_SPEED_M_PER_S = 299792458.0
AU_PER_PC = 206264.80624709636

DEFAULT_ASTRO_BANDS: dict[str, tuple[float, float]] = {
    "FUV": (0.0912, 0.2066),
    "Optical": (0.4, 0.7),
    "Near-IR": (0.7, 5.0),
    "Mid-IR": (5.0, 25.0),
    "Far-IR": (25.0, 500.0),
}

DEFAULT_SLICE_WAVELENGTHS_MICRON: tuple[float, ...] = (0.15, 0.55, 2.0, 10.0, 100.0)

_LENGTH_TO_M = {
    "m": 1.0,
    "meter": 1.0,
    "meters": 1.0,
    "cm": 1e-2,
    "au": 149597870700.0,
    "pc": 149597870700.0 * AU_PER_PC,
    "kpc": 149597870700.0 * AU_PER_PC * 1e3,
}


@dataclass(frozen=True)
class RadiationFieldCube:
    """Container for a SKIRT RadiationFieldProbe FITS cube."""

    data: np.ndarray
    wavelength_micron: np.ndarray
    wavelength_width_micron: np.ndarray | None
    wavelength_left_micron: np.ndarray | None
    wavelength_right_micron: np.ndarray | None
    x: np.ndarray
    y: np.ndarray
    data_unit: str
    spatial_unit: str
    fits_path: Path
    wavelengths_path: Path | None
    header: fits.Header


def load_radiation_field(
    fits_path: str | Path,
    wavelengths_path: str | Path | None = None,
    spatial_unit: str = "pc",
    memmap: bool = True,
) -> RadiationFieldCube:
    """Load a SKIRT RadiationFieldProbe ``*_J.fits`` cube.

    The returned data are ordered as ``(n_lambda, ny, nx)``, matching Astropy's
    FITS convention for a SKIRT cube whose FITS header has axes ``x, y, lambda``.
    Wavelength widths and bin edges are available only when the SKIRT
    ``*_wavelengths.dat`` sidecar file is present.
    """

    fits_path = Path(fits_path)
    if wavelengths_path is None:
        candidate = _infer_wavelength_path(fits_path)
        wavelengths_path = candidate if candidate.exists() else None
    elif wavelengths_path is not None:
        wavelengths_path = Path(wavelengths_path)

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
            raise ValueError(f"Expected a 2D frame or 3D cube in {fits_path}; got shape {data.shape}")

        header = hdul[0].header.copy()
        fits_wavelengths = _read_fits_z_axis(hdul)

    sidecar = _read_wavelength_sidecar(wavelengths_path) if wavelengths_path is not None else None
    if sidecar is not None:
        wavelength, width, left, right = sidecar
    elif fits_wavelengths is not None:
        wavelength = fits_wavelengths
        width = left = right = None
    else:
        raise ValueError(
            f"No wavelength information found for {fits_path}. "
            "Pass wavelengths_path or use a FITS cube with a GRID_POINTS table extension."
        )

    if len(wavelength) != data.shape[0]:
        raise ValueError(
            f"Wavelength axis has {len(wavelength)} bins, but FITS cube has {data.shape[0]} wavelength planes"
        )

    x = _axis_coordinates(header, axis=1, n=data.shape[2], target_unit=spatial_unit)
    y = _axis_coordinates(header, axis=2, n=data.shape[1], target_unit=spatial_unit)

    return RadiationFieldCube(
        data=data,
        wavelength_micron=np.asarray(wavelength, dtype=float),
        wavelength_width_micron=None if width is None else np.asarray(width, dtype=float),
        wavelength_left_micron=None if left is None else np.asarray(left, dtype=float),
        wavelength_right_micron=None if right is None else np.asarray(right, dtype=float),
        x=x,
        y=y,
        data_unit=str(header.get("BUNIT", "")).strip(),
        spatial_unit=spatial_unit,
        fits_path=fits_path,
        wavelengths_path=Path(wavelengths_path) if wavelengths_path is not None else None,
        header=header,
    )


def nearest_wavelength_index(cube: RadiationFieldCube, wavelength_micron: float) -> int:
    """Return the index of the wavelength plane closest to ``wavelength_micron``."""

    if cube.wavelength_micron.size == 0:
        raise ValueError("The cube has no wavelength axis")
    return int(np.nanargmin(np.abs(cube.wavelength_micron - wavelength_micron)))


def to_nu_jnu(cube: RadiationFieldCube) -> np.ndarray:
    """Convert a frequency-style ``J_nu`` cube to ``nu J_nu``.

    S2S currently writes RadiationFieldProbe cubes with SKIRT
    ``fluxOutputStyle="Frequency"``, giving units of ``W/m2/Hz/sr``. This
    function intentionally supports that output style first and raises a clear
    error for unsupported units.
    """

    unit = _normalized_unit(cube.data_unit)
    if "hz" not in unit:
        raise ValueError(
            "to_nu_jnu currently supports frequency-style mean intensity units "
            f"such as 'W/m2/Hz/sr'; got {cube.data_unit!r}"
        )
    wavelength_m = cube.wavelength_micron[:, None, None] * 1e-6
    frequency_hz = LIGHT_SPEED_M_PER_S / wavelength_m
    return np.asarray(cube.data, dtype=float) * frequency_hz


def integrate_band_energy_density(
    cube: RadiationFieldCube,
    lambda_min_micron: float,
    lambda_max_micron: float,
) -> np.ndarray:
    """Integrate radiation energy density over a wavelength interval.

    For frequency-style output, the calculation uses
    ``u = (4*pi/c) * integral J_nu dnu`` and approximates each wavelength bin as
    ``dnu = c * dlambda / lambda**2``.
    """

    width_micron = _require_wavelength_widths(cube)
    unit = _normalized_unit(cube.data_unit)
    if "hz" not in unit:
        raise ValueError(
            "integrate_band_energy_density currently supports frequency-style "
            f"mean intensity units such as 'W/m2/Hz/sr'; got {cube.data_unit!r}"
        )

    mask = _wavelength_mask(cube, lambda_min_micron, lambda_max_micron)
    if not np.any(mask):
        raise ValueError(
            f"No wavelength bins fall in [{lambda_min_micron:g}, {lambda_max_micron:g}] micron"
        )

    wavelength_m = cube.wavelength_micron[mask, None, None] * 1e-6
    width_m = width_micron[mask, None, None] * 1e-6
    jnu = np.asarray(cube.data[mask], dtype=float)
    return np.sum(4.0 * np.pi * jnu * width_m / wavelength_m**2, axis=0)


def plot_wavelength_slices(
    cube: RadiationFieldCube,
    wavelengths_micron: Sequence[float] | None = None,
    quantity: str = "nu_jnu",
    save_path: str | Path | None = None,
    cmap: str = "magma",
    vmin: float | None = None,
    vmax: float | None = None,
):
    """Plot maps at selected wavelengths."""

    if wavelengths_micron is None:
        wavelengths_micron = _default_slice_wavelengths(cube)
    values = _spectral_cube_for_quantity(cube, quantity)
    indices = [nearest_wavelength_index(cube, wavelength) for wavelength in wavelengths_micron]
    frames = [values[index] for index in indices]

    fig, axes = _make_axes(len(frames), figsize_per_panel=(4.2, 3.8))
    norm = _image_norm(frames, vmin=vmin, vmax=vmax, log=True)
    for ax, frame, index in zip(_iter_axes(axes), frames, indices):
        image = _imshow(ax, frame, cube, cmap=cmap, norm=norm)
        ax.set_title(f"{cube.wavelength_micron[index]:.3g} micron")
        _label_spatial_axes(ax, cube)
        fig.colorbar(image, ax=ax, label=_quantity_label(quantity))
    _hide_unused_axes(axes, len(frames))
    _save(fig, save_path)
    return fig, axes


def plot_band_maps(
    cube: RadiationFieldCube,
    bands: Mapping[str, tuple[float, float]] = DEFAULT_ASTRO_BANDS,
    quantity: str = "energy_density",
    save_path: str | Path | None = None,
    cmap: str = "inferno",
    vmin: float | None = None,
    vmax: float | None = None,
):
    """Plot band-integrated radiation-field maps."""

    frames: list[np.ndarray] = []
    labels: list[str] = []
    for label, (lambda_min, lambda_max) in bands.items():
        if quantity != "energy_density":
            raise ValueError("plot_band_maps currently supports quantity='energy_density'")
        frames.append(integrate_band_energy_density(cube, lambda_min, lambda_max))
        labels.append(f"{label} ({lambda_min:g}-{lambda_max:g} micron)")

    fig, axes = _make_axes(len(frames), figsize_per_panel=(4.2, 3.8))
    norm = _image_norm(frames, vmin=vmin, vmax=vmax, log=True)
    for ax, frame, label in zip(_iter_axes(axes), frames, labels):
        image = _imshow(ax, frame, cube, cmap=cmap, norm=norm)
        ax.set_title(label)
        _label_spatial_axes(ax, cube)
        fig.colorbar(image, ax=ax, label=_quantity_label(quantity))
    _hide_unused_axes(axes, len(frames))
    _save(fig, save_path)
    return fig, axes


def plot_spectra(
    cube: RadiationFieldCube,
    apertures: Mapping[str, np.ndarray] | None = None,
    quantity: str = "nu_jnu",
    statistic: str = "median",
    save_path: str | Path | None = None,
):
    """Plot wavelength spectra for full-image or mask-defined apertures."""

    values = _spectral_cube_for_quantity(cube, quantity)
    aperture_masks = _aperture_masks(cube, apertures)

    fig, ax = _plt().subplots(figsize=(6.2, 4.2), constrained_layout=True)
    for label, mask in aperture_masks.items():
        spectrum = _masked_statistic(values, mask, statistic)
        ax.plot(cube.wavelength_micron, spectrum, label=label)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel(_quantity_label(quantity))
    ax.legend()
    _save(fig, save_path)
    return fig, ax


def plot_peak_wavelength_map(
    cube: RadiationFieldCube,
    quantity: str = "nu_jnu",
    save_path: str | Path | None = None,
    cmap: str = "viridis",
):
    """Plot the wavelength where the selected spectral quantity is maximal."""

    values = _spectral_cube_for_quantity(cube, quantity)
    finite_positive = np.isfinite(values) & (values > 0)
    safe_values = np.where(finite_positive, values, -np.inf)
    peak_indices = np.argmax(safe_values, axis=0)
    peak_wavelength = cube.wavelength_micron[peak_indices].astype(float)
    peak_wavelength[~np.any(finite_positive, axis=0)] = np.nan

    fig, ax = _plt().subplots(figsize=(5.2, 4.4), constrained_layout=True)
    image = _imshow(ax, peak_wavelength, cube, cmap=cmap, norm=_colors().LogNorm())
    ax.set_title(f"Peak wavelength of {_quantity_label(quantity)}")
    _label_spatial_axes(ax, cube)
    fig.colorbar(image, ax=ax, label="Peak wavelength [micron]")
    _save(fig, save_path)
    return fig, ax


def plot_view_comparison(
    cubes: Mapping[str, RadiationFieldCube] | Sequence[RadiationFieldCube],
    mode: str = "band",
    band: str | tuple[float, float] = "FUV",
    wavelength_micron: float | None = None,
    quantity: str = "energy_density",
    save_path: str | Path | None = None,
    cmap: str = "inferno",
):
    """Compare the same radiation-field diagnostic across multiple views."""

    cube_items = _cube_items(cubes)
    frames: list[np.ndarray] = []
    titles: list[str] = []
    for label, cube in cube_items:
        if mode == "band":
            lambda_min, lambda_max, band_label = _resolve_band(band)
            if quantity != "energy_density":
                raise ValueError("Band view comparisons currently support quantity='energy_density'")
            frames.append(integrate_band_energy_density(cube, lambda_min, lambda_max))
            titles.append(f"{label}: {band_label}")
        elif mode == "wavelength":
            if wavelength_micron is None:
                raise ValueError("wavelength_micron is required when mode='wavelength'")
            values = _spectral_cube_for_quantity(cube, quantity)
            index = nearest_wavelength_index(cube, wavelength_micron)
            frames.append(values[index])
            titles.append(f"{label}: {cube.wavelength_micron[index]:.3g} micron")
        else:
            raise ValueError("mode must be 'band' or 'wavelength'")

    fig, axes = _make_axes(len(frames), figsize_per_panel=(4.2, 3.8))
    norm = _image_norm(frames, log=True)
    for ax, frame, title, (_, cube) in zip(_iter_axes(axes), frames, titles, cube_items):
        image = _imshow(ax, frame, cube, cmap=cmap, norm=norm)
        ax.set_title(title)
        _label_spatial_axes(ax, cube)
        fig.colorbar(image, ax=ax, label=_quantity_label(quantity))
    _hide_unused_axes(axes, len(frames))
    _save(fig, save_path)
    return fig, axes


def _infer_wavelength_path(fits_path: Path) -> Path:
    name = fits_path.name
    if name.endswith("_J.fits"):
        return fits_path.with_name(name[: -len("_J.fits")] + "_wavelengths.dat")
    return fits_path.with_name(fits_path.stem + "_wavelengths.dat")


def _read_wavelength_sidecar(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    table = np.loadtxt(path, comments="#", ndmin=2)
    if table.shape[1] < 4:
        raise ValueError(f"Expected at least four columns in wavelength table {path}")
    return table[:, 0], table[:, 1], table[:, 2], table[:, 3]


def _read_fits_z_axis(hdul: fits.HDUList) -> np.ndarray | None:
    for hdu in hdul[1:]:
        if hdu.data is None:
            continue
        names = getattr(hdu.data, "names", None) or []
        if names:
            name = "GRID_POINTS" if "GRID_POINTS" in names else names[0]
            return np.asarray(hdu.data[name], dtype=float)
        try:
            return np.asarray(hdu.data.field(0), dtype=float)
        except (AttributeError, IndexError):
            continue
    return None


def _axis_coordinates(header: fits.Header, axis: int, n: int, target_unit: str) -> np.ndarray:
    crpix = float(header.get(f"CRPIX{axis}", (n + 1) / 2.0))
    crval = float(header.get(f"CRVAL{axis}", 0.0))
    cdelt = float(header.get(f"CDELT{axis}", 1.0))
    source_unit = str(header.get(f"CUNIT{axis}", target_unit)).strip()
    coords = (np.arange(n, dtype=float) + 1.0 - crpix) * cdelt + crval
    return _convert_length(coords, source_unit, target_unit)


def _convert_length(values: np.ndarray, source_unit: str, target_unit: str) -> np.ndarray:
    source_key = _normalized_length_unit(source_unit)
    target_key = _normalized_length_unit(target_unit)
    if source_key == target_key:
        return values
    if source_key not in _LENGTH_TO_M or target_key not in _LENGTH_TO_M:
        raise ValueError(f"Cannot convert spatial axis from {source_unit!r} to {target_unit!r}")
    return values * _LENGTH_TO_M[source_key] / _LENGTH_TO_M[target_key]


def _normalized_length_unit(unit: str) -> str:
    return unit.strip().lower()


def _normalized_unit(unit: str) -> str:
    return unit.strip().lower().replace(" ", "")


def _require_wavelength_widths(cube: RadiationFieldCube) -> np.ndarray:
    if cube.wavelength_width_micron is None:
        raise ValueError(
            "Wavelength bin widths are required for band integration. "
            "Load the SKIRT *_wavelengths.dat sidecar file."
        )
    return cube.wavelength_width_micron


def _wavelength_mask(cube: RadiationFieldCube, lambda_min_micron: float, lambda_max_micron: float) -> np.ndarray:
    if lambda_min_micron > lambda_max_micron:
        lambda_min_micron, lambda_max_micron = lambda_max_micron, lambda_min_micron
    return (cube.wavelength_micron >= lambda_min_micron) & (cube.wavelength_micron <= lambda_max_micron)


def _spectral_cube_for_quantity(cube: RadiationFieldCube, quantity: str) -> np.ndarray:
    if quantity == "raw":
        return np.asarray(cube.data, dtype=float)
    if quantity == "nu_jnu":
        return to_nu_jnu(cube)
    raise ValueError("quantity must be 'nu_jnu' or 'raw' for spectral plots")


def _default_slice_wavelengths(cube: RadiationFieldCube) -> list[float]:
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


def _imshow(ax, frame: np.ndarray, cube: RadiationFieldCube, cmap: str, norm=None):
    extent = [cube.x.min(), cube.x.max(), cube.y.min(), cube.y.max()]
    return ax.imshow(frame, origin="lower", extent=extent, aspect="equal", cmap=cmap, norm=norm)


def _image_norm(
    frames: Sequence[np.ndarray],
    vmin: float | None = None,
    vmax: float | None = None,
    log: bool = True,
):
    values = np.concatenate([np.ravel(frame[np.isfinite(frame)]) for frame in frames])
    if values.size == 0:
        return None
    if log:
        positive = values[values > 0]
        if positive.size == 0:
            return _colors().Normalize(vmin=vmin, vmax=vmax)
        return _colors().LogNorm(
            vmin=float(vmin) if vmin is not None else float(np.nanmin(positive)),
            vmax=float(vmax) if vmax is not None else float(np.nanmax(positive)),
        )
    return _colors().Normalize(
        vmin=float(vmin) if vmin is not None else float(np.nanmin(values)),
        vmax=float(vmax) if vmax is not None else float(np.nanmax(values)),
    )


def _label_spatial_axes(ax, cube: RadiationFieldCube) -> None:
    ax.set_xlabel(f"x [{cube.spatial_unit}]")
    ax.set_ylabel(f"y [{cube.spatial_unit}]")


def _quantity_label(quantity: str) -> str:
    if quantity == "nu_jnu":
        return r"$\nu J_\nu$ [W m$^{-2}$ sr$^{-1}$]"
    if quantity == "energy_density":
        return r"Band radiation energy density [J m$^{-3}$]"
    if quantity == "raw":
        return "Mean intensity"
    return quantity


def _save(fig, save_path: str | Path | None) -> None:
    if save_path is None:
        return
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=200, bbox_inches="tight")


def _aperture_masks(
    cube: RadiationFieldCube,
    apertures: Mapping[str, np.ndarray] | None,
) -> dict[str, np.ndarray]:
    if apertures is None:
        return {"All pixels": np.ones(cube.data.shape[1:], dtype=bool)}
    result = {}
    for label, mask in apertures.items():
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != cube.data.shape[1:]:
            raise ValueError(f"Aperture {label!r} has shape {mask.shape}; expected {cube.data.shape[1:]}")
        result[label] = mask
    return result


def _masked_statistic(values: np.ndarray, mask: np.ndarray, statistic: str) -> np.ndarray:
    selected = values[:, mask]
    if selected.size == 0:
        raise ValueError("Aperture mask selects no pixels")
    if statistic == "median":
        return np.nanmedian(selected, axis=1)
    if statistic == "mean":
        return np.nanmean(selected, axis=1)
    if statistic == "sum":
        return np.nansum(selected, axis=1)
    raise ValueError("statistic must be 'median', 'mean', or 'sum'")


def _cube_items(cubes: Mapping[str, RadiationFieldCube] | Sequence[RadiationFieldCube]):
    if isinstance(cubes, Mapping):
        return list(cubes.items())
    return [(f"view {index + 1}", cube) for index, cube in enumerate(cubes)]


def _resolve_band(band: str | tuple[float, float]) -> tuple[float, float, str]:
    if isinstance(band, str):
        if band not in DEFAULT_ASTRO_BANDS:
            raise ValueError(f"Unknown band {band!r}; expected one of {list(DEFAULT_ASTRO_BANDS)}")
        lambda_min, lambda_max = DEFAULT_ASTRO_BANDS[band]
        return lambda_min, lambda_max, band
    lambda_min, lambda_max = band
    return lambda_min, lambda_max, f"{lambda_min:g}-{lambda_max:g} micron"


def _plt():
    import matplotlib.pyplot as plt

    return plt


def _colors():
    from matplotlib import colors

    return colors
