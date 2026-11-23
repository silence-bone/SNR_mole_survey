#!/usr/bin/env python3
"""Make a W44 radio-continuum + H2 2.122 micron RGB overlay.

Radio emission is mapped from black to blue, H2 from black to red, and the
four 30-arcsec telescope regions are drawn in cyan.  The program accepts
either a 2-D radio composite (fast for experimenting with display limits) or
the original 3-D/4-D radio cube, in which case a frequency/channel subset is
median-combined after automatic bad-plane rejection.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import warnings

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / "tmp/mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
from scipy.ndimage import map_coordinates

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.wcs import FITSFixedWarning, WCS


DEFAULT_H2 = "w44_h2_ksub.fits"
DEFAULT_RADIO = "w44_radio_composite.fits"

# Formal J2000 telescope pointings; every circle has a 30 arcsec diameter.
REGIONS = {
    "W44N": SkyCoord("18h55m39.800s", "+01d37m56.400s", frame="icrs"),
    "W44C": SkyCoord("18h56m14.700s", "+01d20m52.700s", frame="icrs"),
    "W44W": SkyCoord("18h55m17.100s", "+01d21m50.000s", frame="icrs"),
    "W44S": SkyCoord("18h56m35.460s", "+01d06m44.100s", frame="icrs"),
}

LABEL_OFFSETS = {
    "W44N": (+93, +41),
    "W44C": (0, +93),
    "W44W": (-95, +55),
    "W44S": (+72, -63),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--h2", default=DEFAULT_H2, help="H2 FITS image")
    parser.add_argument(
        "--radio",
        default=DEFAULT_RADIO,
        help="2-D radio composite or original 3-D/4-D radio cube",
    )
    parser.add_argument("--output", default="w44_radio_h2_overlay.png")
    parser.add_argument(
        "--pdf-output", default="", help="Optional PDF output; empty disables it"
    )
    parser.add_argument(
        "--binning", type=int, default=4, help="Display every Nth H2 pixel"
    )

    frequency = parser.add_argument_group("radio frequency selection")
    frequency.add_argument("--freq-min-mhz", type=float, default=None)
    frequency.add_argument("--freq-max-mhz", type=float, default=None)
    frequency.add_argument(
        "--channels",
        default="",
        help="Channel list/ranges, e.g. '0,2-7,10-15'; overrides frequency limits",
    )
    frequency.add_argument(
        "--list-channels",
        action="store_true",
        help="Print channel centre frequencies and exit",
    )
    frequency.add_argument(
        "--keep-bad-channels",
        action="store_true",
        help="Disable automatic rejection of corrupt/saturated radio planes",
    )

    display = parser.add_argument_group("display scaling")
    display.add_argument("--radio-vmin", type=float, default=None)
    display.add_argument("--radio-vmax", type=float, default=None)
    display.add_argument("--h2-vmin", type=float, default=None)
    display.add_argument("--h2-vmax", type=float, default=None)
    display.add_argument("--radio-percentiles", nargs=2, type=float, default=(5.0, 99.7))
    display.add_argument("--h2-percentiles", nargs=2, type=float, default=(50.0, 99.7))
    display.add_argument(
        "--radio-stretch", choices=("linear", "sqrt", "asinh"), default="asinh"
    )
    display.add_argument(
        "--h2-stretch", choices=("linear", "sqrt", "asinh"), default="asinh"
    )
    display.add_argument(
        "--asinh-softening",
        type=float,
        default=8.0,
        help="Larger values reveal fainter emission in asinh mode",
    )
    display.add_argument("--radio-weight", type=float, default=1.0)
    display.add_argument("--h2-weight", type=float, default=1.0)
    display.add_argument("--radio-color", default="#0066ff")
    display.add_argument("--h2-color", default="#ff2400")
    display.add_argument("--region-color", default="cyan")
    display.add_argument("--colorbars", action="store_true")
    display.add_argument("--title", default="W44 radio continuum + H$_2$ 2.122 $\\mu$m")
    return parser.parse_args()


def decimated_wcs(full_wcs: WCS, step: int) -> WCS:
    """Return WCS appropriate for data sampled as image[::step, ::step]."""
    result = full_wcs.deepcopy()
    result.wcs.crpix = (result.wcs.crpix - 1.0) / step + 1.0
    if result.wcs.has_cd():
        result.wcs.cd = result.wcs.cd * step
    else:
        result.wcs.cdelt = result.wcs.cdelt * step
    return result


def channel_count(data: np.ndarray) -> int:
    if data.ndim == 4:
        return data.shape[1]
    if data.ndim == 3:
        return data.shape[0]
    if data.ndim == 2:
        return 1
    raise ValueError(f"Unsupported radio data shape {data.shape}; expected 2-D, 3-D, or 4-D")


def radio_plane(data: np.ndarray, channel: int) -> np.ndarray:
    if data.ndim == 4:
        return data[0, channel]
    if data.ndim == 3:
        return data[channel]
    if data.ndim == 2 and channel == 0:
        return data
    raise ValueError(f"Cannot get channel {channel} from radio shape {data.shape}")


def channel_frequencies_mhz(header: fits.Header, nchan: int) -> np.ndarray:
    """Read a linear FITS spectral axis and return channel centres in MHz."""
    if nchan == 1:
        return np.array([np.nan])
    spectral_axis = None
    for axis in range(1, int(header.get("NAXIS", 0)) + 1):
        ctype = str(header.get(f"CTYPE{axis}", "")).upper()
        if any(token in ctype for token in ("FREQ", "SPECLN", "VRAD", "VELO")):
            spectral_axis = axis
            break
    if spectral_axis is None:
        raise ValueError("No spectral/frequency axis was found in the radio FITS header")
    crpix = float(header.get(f"CRPIX{spectral_axis}", 1.0))
    crval = float(header[f"CRVAL{spectral_axis}"])
    cdelt = float(header[f"CDELT{spectral_axis}"])
    values = crval + (np.arange(nchan) + 1.0 - crpix) * cdelt
    unit = str(header.get(f"CUNIT{spectral_axis}", "")).strip().lower()
    # This W44 mosaic omits CUNIT3 but stores its ~1.36e9 values in Hz.
    if unit in ("ghz",):
        return values * 1.0e3
    if unit in ("mhz",):
        return values
    if unit in ("khz",):
        return values * 1.0e-3
    if unit in ("hz", "") and np.nanmedian(np.abs(values)) > 1.0e6:
        return values * 1.0e-6
    raise ValueError(f"Cannot interpret spectral unit {unit!r} as frequency")


def parse_channel_spec(spec: str, nchan: int) -> list[int]:
    selected: set[int] = set()
    for item in filter(None, re.split(r"[ ,]+", spec.strip())):
        if "-" in item:
            first, last = map(int, item.split("-", 1))
            selected.update(range(min(first, last), max(first, last) + 1))
        else:
            selected.add(int(item))
    invalid = sorted(index for index in selected if not 0 <= index < nchan)
    if invalid:
        raise ValueError(f"Channel indices outside 0..{nchan - 1}: {invalid}")
    return sorted(selected)


def h2_footprint_in_radio_pixels(
    h2_wcs: WCS, h2_shape: tuple[int, int], radio_wcs: WCS
) -> tuple[np.ndarray, np.ndarray]:
    ny, nx = h2_shape
    samples = 128
    edge_x = np.concatenate(
        [np.linspace(0, nx - 1, samples), np.full(samples, nx - 1),
         np.linspace(nx - 1, 0, samples), np.zeros(samples)]
    )
    edge_y = np.concatenate(
        [np.zeros(samples), np.linspace(0, ny - 1, samples),
         np.full(samples, ny - 1), np.linspace(ny - 1, 0, samples)]
    )
    return radio_wcs.world_to_pixel(h2_wcs.pixel_to_world(edge_x, edge_y))


def select_and_combine_radio(
    data: np.ndarray,
    header: fits.Header,
    h2_wcs: WCS,
    h2_shape: tuple[int, int],
    args: argparse.Namespace,
) -> tuple[np.ndarray, WCS, list[int], list[int], np.ndarray]:
    """Crop, quality-check, and median-combine selected radio channels."""
    nchan = channel_count(data)
    frequencies = channel_frequencies_mhz(header, nchan)
    radio_wcs = WCS(header).celestial
    footprint_x, footprint_y = h2_footprint_in_radio_pixels(h2_wcs, h2_shape, radio_wcs)
    ny, nx = data.shape[-2:]
    x0 = max(0, int(np.floor(np.nanmin(footprint_x))) - 2)
    x1 = min(nx, int(np.ceil(np.nanmax(footprint_x))) + 3)
    y0 = max(0, int(np.floor(np.nanmin(footprint_y))) - 2)
    y1 = min(ny, int(np.ceil(np.nanmax(footprint_y))) + 3)

    if args.channels:
        requested = parse_channel_spec(args.channels, nchan)
    else:
        requested = list(range(nchan))
        if nchan > 1 and args.freq_min_mhz is not None:
            requested = [i for i in requested if frequencies[i] >= args.freq_min_mhz]
        if nchan > 1 and args.freq_max_mhz is not None:
            requested = [i for i in requested if frequencies[i] <= args.freq_max_mhz]
    if not requested:
        raise ValueError("The requested frequency/channel selection contains no radio planes")

    accepted: list[int] = []
    rejected: list[int] = []
    planes: list[np.ndarray] = []
    for channel in requested:
        plane = np.array(radio_plane(data, channel)[y0:y1, x0:x1], dtype=np.float32)
        # This mosaic uses values near +/-3 as invalid/saturated sentinels.
        plane[(~np.isfinite(plane)) | (np.abs(plane) >= 2.9)] = np.nan
        finite_fraction = float(np.isfinite(plane).mean())
        p99_abs = float(np.nanpercentile(np.abs(plane), 99)) if np.isfinite(plane).any() else np.inf
        good = args.keep_bad_channels or nchan == 1 or (finite_fraction > 0.8 and p99_abs < 0.1)
        if good:
            accepted.append(channel)
            planes.append(plane)
        else:
            rejected.append(channel)
    if not planes:
        raise RuntimeError("Every requested radio plane failed the quality check")

    composite = np.nanmedian(np.stack(planes), axis=0).astype(np.float32)
    crop_wcs = radio_wcs.deepcopy()
    crop_wcs.wcs.crpix -= [x0, y0]
    return composite, crop_wcs, accepted, rejected, frequencies


def reproject_to_grid(
    source: np.ndarray, source_wcs: WCS, target_wcs: WCS, target_shape: tuple[int, int]
) -> np.ndarray:
    """Bilinearly sample source onto target without requiring `reproject`."""
    ny, nx = target_shape
    output = np.full((ny, nx), np.nan, dtype=np.float32)
    x = np.arange(nx, dtype=float)
    for y0 in range(0, ny, 256):
        y1 = min(ny, y0 + 256)
        xx, yy = np.meshgrid(x, np.arange(y0, y1, dtype=float))
        sky = target_wcs.pixel_to_world(xx, yy)
        source_x, source_y = source_wcs.world_to_pixel(sky)
        sampled = map_coordinates(
            source,
            [source_y, source_x],
            order=1,
            mode="constant",
            cval=np.nan,
            prefilter=False,
        )
        output[y0:y1] = sampled.astype(np.float32)
    return output


def choose_limits(
    data: np.ndarray, explicit_min: float | None, explicit_max: float | None,
    percentiles: tuple[float, float], name: str
) -> tuple[float, float]:
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        raise RuntimeError(f"{name} has no finite pixels on the common sky grid")
    auto_min, auto_max = np.percentile(finite, percentiles)
    vmin = auto_min if explicit_min is None else explicit_min
    vmax = auto_max if explicit_max is None else explicit_max
    if vmax <= vmin:
        raise ValueError(f"{name} vmax ({vmax:g}) must exceed vmin ({vmin:g})")
    return float(vmin), float(vmax)


def stretched(data: np.ndarray, vmin: float, vmax: float, mode: str, softening: float) -> np.ndarray:
    scaled = np.clip((data - vmin) / (vmax - vmin), 0.0, 1.0)
    if mode == "sqrt":
        scaled = np.sqrt(scaled)
    elif mode == "asinh":
        if softening <= 0:
            raise ValueError("--asinh-softening must be positive")
        scaled = np.arcsinh(softening * scaled) / np.arcsinh(softening)
    return np.nan_to_num(scaled, nan=0.0)


def hex_rgb(value: str) -> np.ndarray:
    value = value.lstrip("#")
    if len(value) != 6:
        raise ValueError(f"Color must be #RRGGBB, got {value!r}")
    return np.array([int(value[i:i + 2], 16) for i in (0, 2, 4)], dtype=float) / 255.0


def add_regions(ax: plt.Axes, color: str) -> None:
    for name, centre in REGIONS.items():
        ax.add_patch(
            SphericalCircle(
                (centre.ra, centre.dec),
                radius=15 * u.arcsec,
                edgecolor=color,
                facecolor="none",
                linewidth=1.8,
                transform=ax.get_transform("icrs"),
                zorder=6,
            )
        )
        dlon, dlat = LABEL_OFFSETS[name]
        position = centre.spherical_offsets_by(dlon * u.arcsec, dlat * u.arcsec)
        label = ax.text(
            position.ra.deg,
            position.dec.deg,
            name,
            transform=ax.get_transform("icrs"),
            color=color,
            fontsize=12,
            fontweight="bold",
            ha="center",
            va="center",
            zorder=7,
        )
        label.set_path_effects(
            [path_effects.Stroke(linewidth=1.5, foreground="black"), path_effects.Normal()]
        )


def main() -> None:
    args = parse_args()
    if args.binning < 1:
        raise ValueError("--binning must be at least 1")
    if args.radio_weight < 0 or args.h2_weight < 0:
        raise ValueError("Image weights cannot be negative")
    warnings.simplefilter("ignore", FITSFixedWarning)

    with fits.open(args.radio, memmap=True) as radio_hdul:
        radio_header = radio_hdul[0].header.copy()
        radio_data = radio_hdul[0].data
        if radio_data is None:
            raise ValueError(f"No image data in {args.radio}")
        nchan = channel_count(radio_data)
        frequencies = channel_frequencies_mhz(radio_header, nchan)
        if args.list_channels:
            if nchan == 1:
                print("The radio input is a 2-D image; it has no selectable frequency channels.")
            else:
                for index, frequency in enumerate(frequencies):
                    print(f"channel {index:2d}: {frequency:10.3f} MHz")
            return

        with fits.open(args.h2, memmap=True) as h2_hdul:
            h2_header = h2_hdul[0].header.copy()
            h2_full = h2_hdul[0].data
            if h2_full is None or h2_full.ndim != 2:
                raise ValueError(f"Expected a 2-D H2 image, got {getattr(h2_full, 'shape', None)}")
            h2_wcs_full = WCS(h2_header).celestial
            h2_image = np.asarray(h2_full[::args.binning, ::args.binning], dtype=np.float32)
            display_wcs = decimated_wcs(h2_wcs_full, args.binning)

            radio_composite, radio_wcs, accepted, rejected, frequencies = select_and_combine_radio(
                radio_data, radio_header, h2_wcs_full, h2_full.shape, args
            )
            radio_image = reproject_to_grid(
                radio_composite, radio_wcs, display_wcs, h2_image.shape
            )

    radio_limits = choose_limits(
        radio_image, args.radio_vmin, args.radio_vmax,
        tuple(args.radio_percentiles), "radio"
    )
    h2_limits = choose_limits(
        h2_image, args.h2_vmin, args.h2_vmax, tuple(args.h2_percentiles), "H2"
    )
    radio_scaled = stretched(
        radio_image, *radio_limits, args.radio_stretch, args.asinh_softening
    )
    h2_scaled = stretched(h2_image, *h2_limits, args.h2_stretch, args.asinh_softening)
    rgb = np.clip(
        args.radio_weight * radio_scaled[..., None] * hex_rgb(args.radio_color)
        + args.h2_weight * h2_scaled[..., None] * hex_rgb(args.h2_color),
        0.0,
        1.0,
    )

    fig = plt.figure(figsize=(9.3, 8.2), dpi=220)
    ax = fig.add_subplot(111, projection=display_wcs)
    ax.imshow(rgb, origin="lower", interpolation="nearest", rasterized=True)
    add_regions(ax, args.region_color)
    ax.set_title(args.title, fontsize=15, pad=11)
    ra_axis, dec_axis = ax.coords[0], ax.coords[1]
    ra_axis.set_axislabel("RA (J2000)", fontsize=13)
    dec_axis.set_axislabel("Dec (J2000)", fontsize=13)
    ra_axis.set_major_formatter("hh:mm:ss")
    dec_axis.set_major_formatter("dd:mm:ss")
    ra_axis.set_ticks(spacing=7.5 * u.arcmin)
    dec_axis.set_ticks(spacing=5 * u.arcmin)
    ra_axis.set_ticklabel(size=10)
    dec_axis.set_ticklabel(size=10, simplify=False)

    if args.colorbars:
        radio_cmap = LinearSegmentedColormap.from_list("radio_black_blue", ["black", args.radio_color])
        h2_cmap = LinearSegmentedColormap.from_list("h2_black_red", ["black", args.h2_color])
        radio_bar = fig.colorbar(
            ScalarMappable(norm=Normalize(*radio_limits), cmap=radio_cmap),
            ax=ax, orientation="horizontal", fraction=0.035, pad=0.08, aspect=35
        )
        radio_bar.set_label(f"Radio intensity ({radio_header.get('BUNIT', 'FITS units')})")
        h2_bar = fig.colorbar(
            ScalarMappable(norm=Normalize(*h2_limits), cmap=h2_cmap),
            ax=ax, orientation="horizontal", fraction=0.035, pad=0.12, aspect=35
        )
        h2_bar.set_label(f"H$_2$ intensity ({h2_header.get('BUNIT', 'FITS units')})")

    fig.tight_layout()
    fig.savefig(args.output, bbox_inches="tight", facecolor="white")
    if args.pdf_output:
        fig.savefig(args.pdf_output, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    if nchan == 1:
        print("Radio input: 2-D composite (no frequency selection applied)")
    else:
        selected_text = ", ".join(
            f"{index} ({frequencies[index]:.1f} MHz)" for index in accepted
        )
        print(f"Accepted radio channels: {selected_text}")
        if rejected:
            print(f"Rejected by quality check: {rejected}")
    print(f"Radio display limits: {radio_limits[0]:.6g} to {radio_limits[1]:.6g}")
    print(f"H2 display limits: {h2_limits[0]:.6g} to {h2_limits[1]:.6g}")
    print(f"Saved {args.output}")
    if args.pdf_output:
        print(f"Saved {args.pdf_output}")


if __name__ == "__main__":
    main()
