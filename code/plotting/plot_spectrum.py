"""Load and plot a one-dimensional FITS spectrum."""

from pathlib import Path
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import c
from astropy.io import fits


# --------------------------- User settings ---------------------------
PROJECT_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_DIR / 'code/other'))
from beam_efficiency import resolve_beam_efficiency, temperature_divisor, iram_header_efficiency
FITS_FILE = PROJECT_DIR / "snr_data/3C391SW_12CO43.fits"
OUTPUT_FILE = PROJECT_DIR / "output/3C391SW_12CO43_spec.png"

# Four velocity-window variables, in km/s.
PLOT_VMIN = 0.0
PLOT_VMAX = 200.0
LINE_VMIN = 40.0
LINE_VMAX = 150.0

BASELINE_ORDER = 1
APPLY_BEAM_EFFICIENCY = True
SHOW_FIGURE = True
# --------------------------------------------------------------------


def load_spectrum(
    fits_file: str | Path,
    apply_beam_efficiency: bool = True,
) -> tuple[np.ndarray, np.ndarray, fits.Header]:
    """Load velocity [km/s] and main-beam temperature [K] from a FITS file."""
    with fits.open(fits_file, memmap=False) as hdul:
        header = hdul[0].header.copy()
        tmb = np.asarray(hdul[0].data, dtype=float).squeeze()

    if tmb.ndim != 1:
        raise ValueError(f"Expected a 1-D spectrum, got shape {tmb.shape}")

    pixel = np.arange(tmb.size) + 1.0  # FITS pixel numbers start at 1.
    frequency_offset = (
        header["CRVAL1"]
        + (pixel - header["CRPIX1"]) * header["CDELT1"]
    )
    velocity = (
        header["VELO-LSR"]
        - c.value * frequency_offset / header["RESTFREQ"]
    ) / 1000.0

    if apply_beam_efficiency:
        beam_efficiency = float(header["BEAMEFF"])
        if iram_header_efficiency(header, fits_file) != 1:
            beam_efficiency = temperature_divisor(header, resolve_beam_efficiency(header)['beff'], fits_file)
        if beam_efficiency <= 0:
            raise ValueError(f"Invalid BEAMEFF: {beam_efficiency}")
        tmb = tmb / beam_efficiency

    return velocity, tmb, header


def title_from_filename(fits_file: str | Path) -> str:
    """Convert ``3C391SW_12CO43.fits`` to ``3C391SW 12CO 4-3``."""
    stem = Path(fits_file).stem
    match = re.fullmatch(r"(.+?)_(\d+CO)(\d)(\d)", stem, flags=re.IGNORECASE)
    if match is None:
        return stem
    source, molecule, upper_j, lower_j = match.groups()
    return f"{source} {molecule.upper()} {upper_j}-{lower_j}"


def plot_spectrum(
    fits_file: str | Path,
    plot_vmin: float,
    plot_vmax: float,
    line_vmin: float,
    line_vmax: float,
    *,
    baseline_order: int = 1,
    apply_beam_efficiency: bool = True,
    ax: plt.Axes | None = None,
    figsize: tuple[float, float] = (6, 4),
    save_dpi: int = 300,
    title: str | None = None,
    output_file: str | Path | None = None,
    show: bool = False,
) -> tuple[plt.Figure, plt.Axes, float]:
    """Plot a FITS spectrum and return ``(figure, axes, rms)``.

    The baseline is fitted within the plot window but outside the line window.
    RMS is calculated from the residuals of those baseline channels.
    """
    plot_low, plot_high = sorted((plot_vmin, plot_vmax))
    line_low, line_high = sorted((line_vmin, line_vmax))
    if not (plot_low <= line_low < line_high <= plot_high):
        raise ValueError("The line window must lie inside the plot window")

    velocity, tmb, header = load_spectrum(
        fits_file, apply_beam_efficiency=apply_beam_efficiency
    )
    finite = np.isfinite(velocity) & np.isfinite(tmb)
    visible = finite & (velocity >= plot_low) & (velocity <= plot_high)
    line = visible & (velocity >= line_low) & (velocity <= line_high)
    baseline_channels = visible & ~line

    if baseline_channels.sum() < baseline_order + 2:
        raise ValueError("Too few channels for the baseline fit")

    baseline_model = np.polynomial.Polynomial.fit(
        velocity[baseline_channels],
        tmb[baseline_channels],
        deg=baseline_order,
    )
    baseline = baseline_model(velocity)
    residual = tmb[baseline_channels] - baseline[baseline_channels]
    rms = float(np.sqrt(np.mean(residual**2)))

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    ax.plot(velocity[visible], tmb[visible], color="black", lw=0.8)
    ax.plot(velocity[visible], baseline[visible], "r--", lw=1.0)
    ax.set_xlim(plot_vmin, plot_vmax)
    ax.set_xlabel(r"$V_{\mathrm{LSR}}$ (km s$^{-1}$)")
    ax.set_ylabel(r"$T_{\mathrm{mb}}$ (K)")
    ax.set_title(title or title_from_filename(fits_file))
    ax.text(
        0.98, 0.95, f"RMS = {rms:.4g} K",
        transform=ax.transAxes, ha="right", va="top",
    )
    fig.tight_layout()

    if output_file is not None:
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_file, dpi=save_dpi, bbox_inches="tight")
    if show:
        plt.show()

    return fig, ax, rms


def main() -> None:
    _, _, rms = plot_spectrum(
        FITS_FILE,
        PLOT_VMIN,
        PLOT_VMAX,
        LINE_VMIN,
        LINE_VMAX,
        baseline_order=BASELINE_ORDER,
        apply_beam_efficiency=APPLY_BEAM_EFFICIENCY,
        output_file=OUTPUT_FILE,
        show=SHOW_FIGURE,
    )
    print(f"Baseline RMS: {rms:.6g} K")
    print(f"Saved plot: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
