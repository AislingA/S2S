"""
instrument_outputs.py

"""
import os

import numpy as np
from astropy.io import fits

def load_instrument_cube(fits_path: str) -> dict:
    """
    Load one SKIRT instrument FITS cube.

    Returns a dictionary with:
    - data
    - header
    - wavelengths
    - unit
    - pixel_scale_arcsec
    """
    if not os.path.exists(fits_path):
        raise FileNotFoundError(f"Cannot find FITS file: {fits_path}.")
    
    with fits.open(fits_path) as hdul:
        data = hdul[0].data
        header = hdul[0].header.copy()

        wavelengths = hdul[1].data["GRID_POINTS"]

    unit = header["BUNIT"]
    pixel_scale_arcsec = abs(header["CDELT1"])

    cube = {
        "data": data,
        "header": header,
        "wavelengths": wavelengths,
        "unit": unit,
        "pixel_scale_arcsec": pixel_scale_arcsec
    }

    return cube

def load_instrument_sed(sed_path: str) -> dict:
    """
    Load one SKIRT instrument SED table.

    Returns a dictionary with:
    - wavelength
    - total
    - transparent
    - primary_direct
    - primary_scattered
    - secondary_direct
    - secondary_scattered
    - secondary_transparent
    """
    if not os.path.exists(sed_path):
        raise FileNotFoundError(f"Cannot find SED file: {sed_path}.")
    
    sed_data = np.loadtxt(sed_path)

    sed = {
        "path": sed_path,
        "filename": os.path.basename(sed_path),
        "wavelength": sed_data[:, 0],
        "total": sed_data[:, 1],
        "transparent": sed_data[:, 2],
        "primary_direct": sed_data[:, 3],
        "primary_scattered": sed_data[:, 4],
        "secondary_direct": sed_data[:, 5],
        "secondary_scattered": sed_data[:, 6],
        "secondary_transparent": sed_data[:, 7]
    }


    return sed

def select_band(cube: dict, band_index: int) -> tuple:
    """
    Select one wavelength plane from the cube.

    Returns:
    - image
    - wavelength
    """
    data = cube["data"]
    wavelengths = cube["wavelengths"]

    if band_index < 0 or band_index >= len(wavelengths):
        raise IndexError(f"Band index {band_index} is outside the cube.")
    
    image = data[band_index, :, :]
    wavelength = wavelengths[band_index]

    return image, wavelength

def get_pixel_scale_pc(cube: dict) -> float:
    """
    Compute the physical pixel scale in parsecs per pixel.

    Parameters
    ----------
    cube: dict
        Cube dictionary from load_instrument_cube().

    Returns
    -------
    float
        Pixel scale in pc/pixel.
    """
    header = cube["header"]

    pixel_scale_arcsec = cube["pixel_scale_arcsec"]
    distance_pc = header["DISTANGD"]

    arcsec_to_radian = np.pi / (180.0 * 3600.0)

    pixel_scale_radian = pixel_scale_arcsec * arcsec_to_radian
    pixel_scale_pc = distance_pc * pixel_scale_radian

    return pixel_scale_pc

def get_fov_pc(cube: dict) -> float:
    """
    Compute the image field of view in parsecs.

    Parameters
    ----------
    cube: dict
        Cube dictionary from load_instrument_cube().

    Returns
    -------
    float
        Field of view across one image side, in pc.
    """
    data = cube["data"]

    pixel_scale_pc = get_pixel_scale_pc(cube)
    num_pixels = data.shape[2]

    fov_pc = num_pixels * pixel_scale_pc

    return fov_pc

def get_image_extent_pc(cube: dict) -> list[float]:
    """
    Build the image extent for plotting axes in pc.

    Parameters
    ----------
    cube: dict
        Cube dictionary from load_instrument_cube().

    Returns
    -------
    list
        Matplotlib extent as [xmin, xmax, ymin, ymax].
    """
    fov_pc = get_fov_pc(cube)
    half_fov_pc = fov_pc / 2.0

    extent = [
        -half_fov_pc,
        half_fov_pc,
        -half_fov_pc,
        half_fov_pc,
    ]

    return extent

def get_view_labels(view: str) -> dict:
    """
    Get plot-axis and line-of-sight labels for a triaxial view.

    Parameters
    ----------
    view: str
        View name: "front", "top", or "side".

    Returns
    -------
    dict
        Labels for x_axis, y_axis, and line_of_sight.
    """
    if view == "front":
        labels = {
            "x_axis": "x",
            "y_axis": "y",
            "line_of_sight": "z",
        }
    elif view == "top":
        labels = {
            "x_axis": "x",
            "y_axis": "z",
            "line_of_sight": "y",
        }
    elif view == "side":
        labels = {
            "x_axis": "y",
            "y_axis": "z",
            "line_of_sight": "x",
        }
    else:
        raise ValueError("view must be 'front', 'top', or 'side'.")

    return labels

def make_metadata_title(
    plot_name: str,
    metadata: dict,
    instrument: str | None = None,
    band_label: str | None = None,
    wavelength_micron: float | None = None,
) -> str:
    """
    Build a consistent plot title from metadata.

    Parameters
    ----------
    plot_name: str
        Main plot description.
    metadata: dict
        Dictionary with time, line-of-sight, and zoom labels.
    instrument: str | None
        Optional instrument label.
    band_label: str | None
        Optional band label.
    wavelength_micron: float | None
        Optional wavelength in micron.

    Returns
    -------
    str
        Formatted plot title.
    """
    title_parts = [plot_name]

    if instrument is not None:
        title_parts.append(instrument)

    if band_label is not None:
        title_parts.append(band_label)

    if wavelength_micron is not None:
        title_parts.append(f"{wavelength_micron:.3f} micron")

    main_title = " ".join(title_parts)

    detail_parts = []

    if "snapshot_time_myr" in metadata:
        detail_parts.append(f"t = {metadata['snapshot_time_myr']:.2f} Myr")

    if "line_of_sight" in metadata:
        detail_parts.append(f"LOS = {metadata['line_of_sight']}")

    if "zoom_label" in metadata:
        detail_parts.append(metadata["zoom_label"])

    if detail_parts:
        details = ", ".join(detail_parts)
        return f"{main_title}: {details}"

    return main_title

def apply_gaussian_psf(rgb_cube: np.ndarray, fwhm_pixels: float) -> np.ndarray:
    """
    Apply a simple Gaussian PSF to each RGB channel.

    Parameters
    ----------
    rgb_cube: np.ndarray
        RGB cube with shape (3, ny, nx).
    fwhm_pixels: float
        Gaussian full width at half maximum in pixels.

    Returns
    -------
    np.ndarray
        Smoothed RGB cube with shape (3, ny, nx).
    """
    from astropy.convolution import Gaussian2DKernel, convolve

    if fwhm_pixels <= 0:
        return rgb_cube

    fwhm_to_sigma = 2.0 * np.sqrt(2.0 * np.log(2.0))
    sigma_pixels = fwhm_pixels / fwhm_to_sigma

    kernel = Gaussian2DKernel(x_stddev=sigma_pixels)

    smoothed_channels = []

    for channel in rgb_cube:
        smoothed_channel = convolve(
            channel,
            kernel,
            normalize_kernel=True,
            boundary="extend"
        )
        smoothed_channels.append(smoothed_channel)

    smoothed_cube = np.stack(smoothed_channels)

    return smoothed_cube

def get_rgb_preset(preset_name: str) -> dict:
    """
    Get RGB band choices and labels for a named RGB preset.

    Parameters
    ----------
    preset_name: str
        Either "observation" or "warm_dust".

    Returns
    -------
    dict
        RGB indices, band labels, and title text.
    """
    if preset_name == "observation":
        preset = {
            "red_index": 6,
            "green_index": 2,
            "blue_index": 0,
            "red_label": "F444W",
            "green_label": "F200W",
            "blue_label": "F090W",
            "title": "JWST observation RGB: R=F444W, G=F200W, B=F090W",
        }
    elif preset_name == "warm_dust":
        preset = {
            "red_index": 8,
            "green_index": 7,
            "blue_index": 6,
            "red_label": "F1500W",
            "green_label": "F770W",
            "blue_label": "F444W",
            "title": "JWST warm-dust RGB: R=F1500W, G=F770W, B=F444W",
        }
    else:
        raise ValueError("preset_name must be 'observation' or 'warm_dust'.")

    return preset

def make_rgb_cube(cube: dict, red_index: int = 6, green_index: int = 2, blue_index: int = 0) -> np.ndarray:
    """
    Select three wavelength planes for an RGB image.

    Default JWST bands:
    - red: F444W
    - green: F200W
    - blue: F090W

    Returns
    -------
    np.ndarray
        RGB cube with shape (3, ny, nx).
    """
    red_image, _ = select_band(cube, red_index)
    green_image, _ = select_band(cube, green_index)
    blue_image, _ = select_band(cube, blue_index)

    rgb_cube = np.stack([
        red_image,
        green_image,
        blue_image
    ])

    return rgb_cube

def make_rgb_from_preset(cube: dict, preset_name: str) -> tuple[np.ndarray, dict]:
    """
    Build an RGB cube using a named RGB preset.

    Parameters
    ----------
    cube: dict
        Cube dictionary from load_instrument_cube().
    preset_name: str
        Either "observation" or "warm_dust".

    Returns
    -------
    tuple
        rgb_cube, preset
    """
    preset = get_rgb_preset(preset_name)

    rgb_cube = make_rgb_cube(
        cube,
        red_index=preset["red_index"],
        green_index=preset["green_index"],
        blue_index=preset["blue_index"],
    )

    return rgb_cube, preset

def stretch_rgb(rgb_cube: np.ndarray, lower_percent: float = 1.0, upper_percent: float = 99.0) -> np.ndarray:
    """
    Stretch an RGB cube so it can be displayed as an image.

    Each channel is stretched separately.

    Parameters
    ----------
    rgb_cube: np.ndarray
        RGB cube with shape (3, ny, nx).
    lower_percent: float
        Low percentile used as the dark value.
    upper_percent: float
        High percentile used as the bright value.

    Returns
    -------
    np.ndarray
        RGB image with shape (ny, nx, 3), with values between 0 and 1.
    """
    stretched_channels = []

    for channel in rgb_cube:
        positive_values = channel[channel > 0]

        if len(positive_values) == 0:
            stretched_channel = np.zeros_like(channel)
        else:
            vmin = np.percentile(positive_values, lower_percent)
            vmax = np.percentile(positive_values, upper_percent)

            stretched_channel = (channel - vmin) / (vmax - vmin)
            stretched_channel = np.clip(stretched_channel, 0, 1)

        stretched_channels.append(stretched_channel)

    rgb_image = np.dstack(stretched_channels)

    return rgb_image

def stretch_rgb_asinh(
    rgb_cube: np.ndarray,
    lower_percent: float = 1.0,
    upper_percent: float = 98.0,
    asinh_a: float = 0.01,
) -> np.ndarray:
    """
    Stretch an RGB cube with an asinh stretch.

    This is useful for astronomy images because it reveals faint structure while
    allowing the brightest compact sources to saturate.

    Parameters
    ----------
    rgb_cube: np.ndarray
        RGB cube with shape (3, ny, nx).
    lower_percent: float
        Low percentile used as the dark value.
    upper_percent: float
        High percentile used as the bright value.
    asinh_a: float
        Controls how aggressive the asinh stretch is. Smaller values are more aggressive.

    Returns
    -------
    np.ndarray
        RGB image with shape (ny, nx, 3).
    """
    from astropy.visualization import AsymmetricPercentileInterval, ImageNormalize, AsinhStretch

    stretched_channels = []

    for channel in rgb_cube:
        positive_values = channel[channel > 0]

        if len(positive_values) == 0:
            stretched_channel = np.zeros_like(channel)
        else:
            interval = AsymmetricPercentileInterval(lower_percent, upper_percent)
            vmin, vmax = interval.get_limits(positive_values)

            norm = ImageNormalize(
                vmin=vmin,
                vmax=vmax,
                stretch=AsinhStretch(a=asinh_a),
            )

            stretched_channel = norm(channel)
            stretched_channel = np.ma.filled(stretched_channel, 0)
            stretched_channel = np.clip(stretched_channel, 0, 1)

        stretched_channels.append(stretched_channel)

    rgb_image = np.dstack(stretched_channels)

    return rgb_image

def plot_rgb_image(
    rgb_image: np.ndarray,
    save_path: str | None = None,
    title: str | None = None,
    snapshot_name: str | None = None,
    snapshot_time_myr: float | None = None,
    distance_pc: float | None = None,
    fov_pc: float | None = None,
    view: str | None = None,
    psf_label: str | None = None,
    stretch_label: str | None = None
):
    """
    Plot a stretched RGB image with optional labels.

    Parameters
    ----------
    rgb_image: np.ndarray
        RGB image with shape (ny, nx, 3).
    save_path: str | None
        Optional path for saving the figure.
    title: str | None
        Optional title shown above the image.
    snapshot_name: str | None
        Optional snapshot label.
    snapshot_time_myr: float | None
        Optional snapshot time in Myr.
    distance_pc: float | None
        Optional instrument distance in pc.
    fov_pc: float | None
        Optional field of view in pc.
    view: str | None
        Optional 

    Returns
    -------
    tuple
        fig, ax
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 6))

    ax.imshow(rgb_image, origin="lower")
    ax.axis("off")

    if title is not None:
        ax.set_title(title)

    info_lines = []

    if snapshot_name is not None:
        info_lines.append(snapshot_name)

    if snapshot_time_myr is not None:
        info_lines.append(f"t = {snapshot_time_myr:.2f} Myr")

    if distance_pc is not None:
        info_lines.append(f"distance = {distance_pc:.0f} pc")

    if fov_pc is not None:
        info_lines.append(f"FOV = {fov_pc:.2f} pc")

    if view is not None:
        labels = get_view_labels(view)
        info_lines.append(f"LOS = {labels['line_of_sight']}")

    if psf_label is not None:
        info_lines.append(psf_label)

    if stretch_label is not None:
        info_lines.append(stretch_label)

    if info_lines:
        info_text = "\n".join(info_lines)

        ax.text(
            0.03,
            0.97,
            info_text,
            color="white",
            fontsize=9,
            ha="left",
            va="top",
            transform=ax.transAxes,
            bbox={
                "facecolor": "black",
                "alpha": 0.6,
                "edgecolor": "none",
            },
        )

    if save_path is not None:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")

    return fig, ax

def add_view_arrows(ax, view: str) -> None:
    """
    Add simple triaxial orientation labels to a plot.

    Parameters
    ----------
    ax:
        Matplotlib axis to draw on.
    view: str
        View name: "front", "top", or "side".
    """
    labels = get_view_labels(view)

    ax.annotate(
        labels["x_axis"],
        xy=(0.18, 0.08),
        xytext=(0.08, 0.08),
        xycoords="axes fraction",
        textcoords="axes fraction",
        color="white",
        arrowprops={
            "arrowstyle": "->",
            "color": "white",
        },
    )

    ax.annotate(
        labels["y_axis"],
        xy=(0.08, 0.18),
        xytext=(0.08, 0.08),
        xycoords="axes fraction",
        textcoords="axes fraction",
        color="white",
        arrowprops={
            "arrowstyle": "->",
            "color": "white",
        },
    )

    ax.text(
        0.03,
        0.03,
        f"LOS: {labels['line_of_sight']}",
        color="white",
        fontsize=9,
        ha="left",
        va="bottom",
        transform=ax.transAxes,
    )

def plot_sed(sed: dict, save_path: str | None = None):
    """
    Plot the total SED.

    Parameters
    ----------
    sed: dict
        SED dictionary from load_instrument_sed().
    save_path: str | None
        Optional path for saving the figure.

    Returns
    -------
    tuple
        fig, ax
    """
    import matplotlib.pyplot as plt

    wavelength = sed["wavelength"]
    total_flux = sed["total"]

    fig, ax = plt.subplots(figsize=(6,4))

    ax.plot(wavelength, total_flux, marker="o")
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel("Flux Density [Jy]")
    filename = sed.get("filename", "unknown SED file")
    ax.set_title(f"Total SED: {filename}")


    if save_path is not None:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")

    return fig, ax

def plot_component_seds(sed: dict, save_path: str | None = None):
    """
    Plot the total SED and the main component SEDs.

    Parameters
    ----------
    sed: dict
        SED dictionary from load_instrument_sed().
    save_path: str | None
        Optional path for saving the figure.

    Returns
    -------
    tuple
        fig, ax
    """
    import matplotlib.pyplot as plt

    wavelength = sed["wavelength"]

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.plot(wavelength, sed["total"], marker="o", label="total")
    ax.plot(wavelength, sed["primary_direct"], marker="o", label="primary direct")
    ax.plot(wavelength, sed["primary_scattered"], marker="o", label="primary scattered")
    ax.plot(wavelength, sed["secondary_direct"], marker="o", label="secondary direct")
    ax.plot(wavelength, sed["secondary_scattered"], marker="o", label="secondary scattered")

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("Wavelength [micron]")
    ax.set_ylabel("Flux density [Jy]")
    filename = sed.get("filename", "unknown SED file")
    ax.set_title(f"Component SEDs: {filename}")
    ax.legend()

    if save_path is not None:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")

    return fig, ax

def plot_band_image(cube: dict, band_index: int, save_path: str | None = None, view: str = "front"):
    """
    Plot one wavelength plane from an instrument cube.

    Parameters
    ----------
    cube: dict
        Cube dictionary from load_instrument_cube().
    band_index: int
        Index of the wavelength plane to plot.
    save_path: str | None
        Optional path for saving the figure.

    Returns
    -------
    tuple
        fig, ax
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    image, wavelength = select_band(cube, band_index)

    extent = get_image_extent_pc(cube)
    labels = get_view_labels(view)

    positive_values = image[image > 0]
    vmin = np.percentile(positive_values, 1)
    vmax = np.percentile(positive_values, 99)

    fig, ax = plt.subplots(figsize=(6, 5))

    plot = ax.imshow(
        image,
        origin="lower",
        extent=extent,
        norm=LogNorm(vmin=vmin, vmax=vmax),
        cmap="magma"
    )

    ax.set_title(f"Band {band_index}: {wavelength:.3f} micron")
    ax.set_xlabel(f"{labels['x_axis']} [pc]")
    ax.set_ylabel(f"{labels['y_axis']} [pc]")
    ax.text(
        0.03,
        0.03,
        f"LOS: {labels['line_of_sight']}",
        color="white",
        fontsize=9,
        ha="left",
        va="bottom",
        transform=ax.transAxes,
    )

    fig.colorbar(plot, ax=ax, label=cube["unit"])

    if save_path is not None:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")

    return fig, ax

def load_component_cubes(base_path: str) -> dict:
    """
    Load the main SKIRT instrument component cubes.

    Parameters
    ----------
    base_path: str
        File path without the component name and .fits ending.
        Example: "outputs/snapshot_150_jwst_front"

    Returns
    -------
    dict
        Dictionary of loaded component cubes.
    """
    component_names = [
        "total",
        "transparent",
        "primarydirect",
        "primaryscattered",
        "secondarydirect",
        "secondaryscattered",
    ]

    cubes = {}

    for component_name in component_names:
        fits_path = f"{base_path}_{component_name}.fits"
        cubes[component_name] = load_instrument_cube(fits_path)

    return cubes

def plot_component_maps(component_cubes: dict, band_index: int, save_path: str | None = None):
    """
    Plot several component maps for one wavelength band.

    Parameters
    ----------
    component_cubes: dict
        Dictionary from load_component_cubes().
    band_index: int
        Index of the wavelength plane to plot.
    save_path: str | None
        Optional path for saving the figure.

    Returns
    -------
    tuple
        fig, axes
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    component_names = [
        "total",
        "transparent",
        "primarydirect",
        "primaryscattered",
        "secondarydirect",
        "secondaryscattered",
    ]

    fig, axes = plt.subplots(2, 3, figsize=(12, 8))

    for ax, component_name in zip(axes.ravel(), component_names):
        cube = component_cubes[component_name]
        image, wavelength = select_band(cube, band_index)

        positive_values = image[image > 0]
        vmin = np.percentile(positive_values, 1)
        vmax = np.percentile(positive_values, 99)

        plot = ax.imshow(
            image,
            origin="lower",
            norm=LogNorm(vmin=vmin, vmax=vmax),
            cmap="magma"
        )

        ax.set_title(component_name)
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(f"Component maps at {wavelength:.3f} micron")
    fig.colorbar(plot, ax=axes.ravel().tolist(), label=component_cubes["total"]["unit"])

    if save_path is not None:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")

    return fig, axes

def make_ratio_map(numerator: np.ndarray, denominator: np.ndarray, min_denominator: float = 1e-30) -> np.ndarray:
    """
    Safely divide two images to make a ratio map.

    Parameters
    ----------
    numerator: np.ndarray
        Image used on the top of the ratio.
    denominator: np.ndarray
        Image used on the bottom of the ratio.
    min_denominator: float
        Denominator values smaller than this are masked.

    Returns
    -------
    np.ndarray
        Ratio image. Bad pixels are set to np.nan.
    """
    ratio = np.full_like(numerator, np.nan, dtype=float)

    good_pixels = denominator > min_denominator

    ratio[good_pixels] = numerator[good_pixels] / denominator[good_pixels]

    return ratio

def plot_ratio_map(
    component_cubes: dict,
    numerator_name: str,
    denominator_name: str,
    band_index: int,
    save_path: str | None = None
):
    """
    Plot a ratio map between two component cubes.

    Parameters
    ----------
    component_cubes: dict
        Dictionary from load_component_cubes().
    numerator_name: str
        Name of the component used on top of the ratio.
    denominator_name: str
        Name of the component used on bottom of the ratio.
    band_index: int
        Index of the wavelength plane to plot.
    save_path: str | None
        Optional path for saving the figure.

    Returns
    -------
    tuple
        fig, ax
    """
    import matplotlib.pyplot as plt

    numerator_cube = component_cubes[numerator_name]
    denominator_cube = component_cubes[denominator_name]

    numerator_image, wavelength = select_band(numerator_cube, band_index)
    denominator_image, _ = select_band(denominator_cube, band_index)

    ratio = make_ratio_map(numerator_image, denominator_image)

    finite_values = ratio[np.isfinite(ratio)]
    vmin = np.percentile(finite_values, 1)
    vmax = np.percentile(finite_values, 99)

    fig, ax = plt.subplots(figsize=(6, 5))

    plot = ax.imshow(
        ratio,
        origin="lower",
        vmin=vmin,
        vmax=vmax,
        cmap="viridis"
    )

    ax.set_title(f"{numerator_name} / {denominator_name} at {wavelength:.3f} micron")
    ax.set_xlabel("x pixel")
    ax.set_ylabel("y pixel")

    fig.colorbar(plot, ax=ax, label="ratio")

    if save_path is not None:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")

    return fig, ax

def image_flux_jy(cube: dict, band_index: int) -> float:
    """
    Integrate one FITS image plane into total flux density.

    This assumes the image unit is MJy/sr.

    Parameters
    ----------
    cube: dict
        Cube dictionary from load_instrument_cube().
    band_index: int
        Index of the wavelength plane to integrate.

    Returns
    -------
    float
        Image-integrated flux density in Jy.
    """
    image, _ = select_band(cube, band_index)

    pixel_scale_arcsec = cube["pixel_scale_arcsec"]
    arcsec_to_radian = np.pi / (180.0 * 3600.0)

    pixel_scale_radian = pixel_scale_arcsec * arcsec_to_radian
    pixel_solid_angle = pixel_scale_radian**2

    flux_jy = np.sum(image) * 1e6 * pixel_solid_angle

    return flux_jy

def compare_image_flux_to_sed(cube: dict, sed: dict) -> dict:
    """
    Compare image-integrated flux with SED flux for every band.

    Parameters
    ----------
    cube: dict
        Cube dictionary from load_instrument_cube().
    sed: dict
        SED dictionary from load_instrument_sed().

    Returns
    -------
    dict
        Comparison arrays for wavelength, image flux, SED flux, and ratio.
    """
    wavelengths = cube["wavelengths"]

    image_fluxes = []

    for band_index in range(len(wavelengths)):
        flux = image_flux_jy(cube, band_index)
        image_fluxes.append(flux)

    image_fluxes = np.array(image_fluxes)
    sed_fluxes = sed["total"]

    ratio = image_fluxes / sed_fluxes

    comparison = {
        "wavelength": wavelengths,
        "image_flux": image_fluxes,
        "sed_flux": sed_fluxes,
        "ratio": ratio,
    }

    return comparison

def main() -> None:
    """
    Run the first-pass JWST front-view diagnostics for snapshot_150.
    """
    snapshot_name = "snapshot_150"
    snapshot_time_myr = 0.12
    view = "front"

    base_path = "outputs/snapshot_150_jwst_front"
    total_fits_path = f"{base_path}_total.fits"
    sed_path = f"{base_path}_sed.dat"

    cube = load_instrument_cube(total_fits_path)
    sed = load_instrument_sed(sed_path)
    component_cubes = load_component_cubes(base_path)

    fwhm_pixels = 1.5
    stretch_label = "stretch = asinh, 1-98%, a=0.01"
    psf_label = f"PSF = Gaussian FWHM {fwhm_pixels:.1f} px"

    for preset_name in ["observation", "warm_dust"]:
        rgb_cube, preset = make_rgb_from_preset(cube, preset_name)
        rgb_cube = apply_gaussian_psf(rgb_cube, fwhm_pixels=fwhm_pixels)
        rgb_image = stretch_rgb_asinh(rgb_cube)

        plot_rgb_image(
            rgb_image,
            save_path=f"outputs/snapshot_150_jwst_front_rgb_{preset_name}.png",
            title=preset["title"],
            snapshot_name=snapshot_name,
            snapshot_time_myr=snapshot_time_myr,
            distance_pc=cube["header"]["DISTANGD"],
            fov_pc=get_fov_pc(cube),
            view=view,
            psf_label=psf_label,
            stretch_label=stretch_label
        )

    plot_sed(
        sed,
        save_path="outputs/snapshot_150_jwst_front_sed.png",
    )

    plot_component_seds(
        sed,
        save_path="outputs/snapshot_150_jwst_front_component_seds.png",
    )

    plot_band_image(
        cube,
        band_index=6,
        save_path="outputs/snapshot_150_jwst_front_band6.png",
        view=view,
    )

    plot_component_maps(
        component_cubes,
        band_index=6,
        save_path="outputs/snapshot_150_jwst_front_component_maps_band6.png",
    )

    plot_ratio_map(
        component_cubes,
        numerator_name="total",
        denominator_name="transparent",
        band_index=6,
        save_path="outputs/snapshot_150_jwst_front_total_over_transparent_band6.png",
    )

    comparison = compare_image_flux_to_sed(cube, sed)

    print("Finished instrument-output diagnostics.")
    print("Image/SED flux ratios:")
    print(comparison["ratio"])

if __name__ == "__main__":
    main()