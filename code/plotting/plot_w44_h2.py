#!/usr/bin/env python3
"""Plot the W44 H2 image with four 30-arcsec observing regions."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import warnings

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[2] / "tmp/mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.wcs import FITSFixedWarning, WCS


DEFAULT_H2 = "w44_h2_ksub.fits"

# Formal J2000 telescope pointing coordinates supplied by the user.
# Every circle below has a fixed diameter of 30 arcsec.
REGIONS = {
    "W44N": SkyCoord("18h55m39.800s", "+01d37m56.400s", frame="icrs"),
    "W44C": SkyCoord("18h56m14.700s", "+01d20m52.700s", frame="icrs"),
    "W44W": SkyCoord("18h55m17.100s", "+01d21m50.000s", frame="icrs"),
    "W44S": SkyCoord("18h56m35.460s", "+01d06m44.100s", frame="icrs"),
}

# Spherical longitude/latitude offsets in arcsec for uncluttered labels.
LABEL_OFFSETS = {
    "W44N": (+93, +41),
    "W44C": (0, +93),
    "W44W": (-95, +55),
    "W44S": (+72, -63),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h2", default=DEFAULT_H2, help="Input H2 FITS image")
    parser.add_argument("--output", default="w44_h2_regions.png")
    parser.add_argument(
        "--pdf-output",
        default="w44_h2_regions.pdf",
        help="PDF output; pass an empty string to disable",
    )
    parser.add_argument(
        "--binning",
        type=int,
        default=4,
        help="Display every Nth H2 pixel to reduce memory use (default: 4)",
    )
    parser.add_argument("--vmin", type=float, default=0.0)
    parser.add_argument("--vmax", type=float, default=1.0e-16)
    return parser.parse_args()


def decimated_wcs(full_wcs: WCS, step: int) -> WCS:
    """Return a WCS for data sampled as image[::step, ::step].

    WCS.slice currently warns and leaves the scale unchanged for this FITS
    header because it uses a CD matrix. Applying the pixel-centre relation
    explicitly avoids that failure mode.
    """
    result = full_wcs.deepcopy()
    result.wcs.crpix = (result.wcs.crpix - 1.0) / step + 1.0
    if result.wcs.has_cd():
        result.wcs.cd = result.wcs.cd * step
    else:
        result.wcs.cdelt = result.wcs.cdelt * step
    return result


def main() -> None:
    args = parse_args()
    if args.binning < 1:
        raise ValueError("--binning must be at least 1")
    if args.vmax <= args.vmin:
        raise ValueError("--vmax must be greater than --vmin")

    warnings.simplefilter("ignore", FITSFixedWarning)
    with fits.open(args.h2, memmap=True) as hdul:
        header = hdul[0].header.copy()
        full_data = hdul[0].data
        if full_data.ndim != 2:
            raise ValueError(f"Expected a 2D H2 image, got shape {full_data.shape}")
        image_data = np.asarray(
            full_data[:: args.binning, :: args.binning], dtype=np.float32
        )
        full_wcs = WCS(header).celestial
        image_wcs = decimated_wcs(full_wcs, args.binning)

    fig = plt.figure(figsize=(9.3, 8.1), dpi=220)
    ax = fig.add_subplot(111, projection=image_wcs)
    image = ax.imshow(
        image_data,
        origin="lower",
        cmap="viridis",
        vmin=args.vmin,
        vmax=args.vmax,
        interpolation="nearest",
        rasterized=True,
    )

    region_color = "#ff4a4a"
    for name, centre in REGIONS.items():
        circle = SphericalCircle(
            (centre.ra, centre.dec),
            radius=15 * u.arcsec,
            edgecolor=region_color,
            facecolor="none",
            linewidth=1.8,
            transform=ax.get_transform("icrs"),
            zorder=6,
        )
        ax.add_patch(circle)
        dlon, dlat = LABEL_OFFSETS[name]
        label_coord = centre.spherical_offsets_by(dlon * u.arcsec, dlat * u.arcsec)
        label = ax.text(
            label_coord.ra.deg,
            label_coord.dec.deg,
            name,
            transform=ax.get_transform("icrs"),
            color=region_color,
            fontsize=13,
            fontweight="bold",
            ha="center",
            va="center",
            zorder=7,
        )
        label.set_path_effects(
            [path_effects.Stroke(linewidth=1.4, foreground="black", alpha=0.25), path_effects.Normal()]
        )

    ra_axis = ax.coords[0]
    dec_axis = ax.coords[1]
    ra_axis.set_axislabel("RA (J2000)", fontsize=14)
    dec_axis.set_axislabel("DEC (J2000)", fontsize=14)
    ra_axis.set_major_formatter("hh:mm:ss")
    dec_axis.set_major_formatter("dd:mm:ss")
    # 7.5 arcmin on the sky corresponds to 30 seconds of RA near the equator.
    ra_axis.set_ticks(spacing=7.5 * u.arcmin)
    dec_axis.set_ticks(spacing=5 * u.arcmin)
    ra_axis.set_ticklabel(size=11)
    dec_axis.set_ticklabel(size=11, simplify=False)

    fig.tight_layout()
    fig.savefig(args.output, bbox_inches="tight")
    if args.pdf_output:
        fig.savefig(args.pdf_output, bbox_inches="tight")
    plt.close(fig)

    finite = image_data[np.isfinite(image_data)]
    print(f"H2 finite range: {finite.min():.6e} to {finite.max():.6e}")
    print(f"H2 median: {np.median(finite):.6e}")
    print(f"Saved figure to {args.output}")
    if args.pdf_output:
        print(f"Saved vector figure to {args.pdf_output}")


if __name__ == "__main__":
    main()
