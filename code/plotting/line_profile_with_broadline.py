from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# ------------------------------------------------------------------
# Configuration: three 12CO lines and the two HCO+ lines currently
# present in this directory. To add a third HCO+ transition later,
# append one dictionary with its key, label, scale, offset and color.
# ------------------------------------------------------------------
broad_source = "w44w"
broad_sigma_min = 5.0        # km/s; strict selection: sigma > 5
broad_vmin, broad_vmax = -50.0, 150.0
broad_root = Path("/Users/silence/Desktop/W44")  # directory containing the FITS and fit txt files

broad_line_settings = [
    dict(key="12CO43", label=r"$^{12}\mathrm{CO}$ 4-3", beff=0.61, scale=1.0,  offset=15.0, color="#D55E00"),
    dict(key="12CO21", label=r"$^{12}\mathrm{CO}$ 2-1", beff=0.59, scale=1.0,  offset=10.0, color="#0072B2"),
    dict(key="12CO10", label=r"$^{12}\mathrm{CO}$ 1-0", beff=0.75, scale=1.0,  offset=5.0,  color="#CC79A7"),
    dict(key="HCOp10", label=r"$\mathrm{HCO}^{+}$ 1-0 $\times 10$", beff=0.78, scale=10.0, offset=0.0,  color="#009E73"),
    dict(key="HCOp32", label=r"$\mathrm{HCO}^{+}$ 3-2 $\times 10$", beff=0.82, scale=10.0, offset=-5.0, color="#E69F00"),
]


def read_class_fits_spectrum(filename, beff=1.0):
    """Read the 1-D CLASS FITS files used here and return v_LSR and T_mb."""
    header = fits.getheader(filename)
    spectrum = np.asarray(fits.getdata(filename), dtype=float).squeeze()
    if spectrum.ndim != 1:
        raise ValueError(f"Expected a 1-D spectrum after squeeze: {filename}")

    pixel = np.arange(spectrum.size, dtype=float) + 1.0  # FITS pixels are 1-based
    if header.get("DELTAV") is not None and header.get("VELO-LSR") is not None:
        velocity = (
            header["VELO-LSR"]
            + (pixel - header["CRPIX1"]) * header["DELTAV"]
        ) / 1e3
    else:
        # Radio definition, used only if the CLASS velocity keywords are absent.
        frequency = header["CRVAL1"] + (pixel - header["CRPIX1"]) * header["CDELT1"]
        rest_frequency = header.get("RESTFRQ", header.get("RESTFREQ"))
        velocity = 299792.458 * (1.0 - frequency / rest_frequency)

    # These FITS spectra have not been beam-efficiency corrected.
    # Convert every observed spectrum to T_mb using T_mb = data / eta_mb.
    spectrum = spectrum / beff
    return velocity, spectrum


def read_broad_fit_components(filename, sigma_min):
    """Return only fit rows whose sigma is strictly larger than sigma_min."""
    values = np.loadtxt(filename, comments="#", ndmin=2)
    components = []
    for row in values:
        if np.isfinite(row[[2, 4, 6]]).all() and row[4] > sigma_min:
            components.append({
                "component": int(row[0]), "type": int(row[1]),
                "v0": float(row[2]), "sigma": float(row[4]),
                "amp": float(row[6]),
            })
    return components


def broad_gaussian(velocity, amp, v0, sigma):
    return amp * np.exp(-0.5 * ((velocity - v0) / sigma) ** 2)


def plot_stacked_broad_components(
    source=broad_source,
    settings=broad_line_settings,
    root=broad_root,
    sigma_min=broad_sigma_min,
    vmin=broad_vmin,
    vmax=broad_vmax,
    outfile=None,
):
    """Stack spectra and overlay only individual Gaussian components above sigma_min."""
    root = Path(root)
    if outfile is None:
        outfile = root / f"{source}_broad_components.png"

    fig, ax = plt.subplots(figsize=(6.0, 8.0), dpi=200)
    selected_rows = []

    for item in settings:
        key = item["key"]
        datafile = root / f"{source}_{key}_re.fits"
        fitfile = root / f"{source}_{key}_fit.txt"
        if not datafile.exists() or not fitfile.exists():
            print(f"[Skip] missing data or fit file for {key}")
            continue

        velocity, spectrum = read_class_fits_spectrum(datafile, item["beff"])
        keep = (
            np.isfinite(velocity) & np.isfinite(spectrum)
            & (velocity >= vmin) & (velocity <= vmax)
        )
        velocity, spectrum = velocity[keep], spectrum[keep]
        order = np.argsort(velocity)
        velocity, spectrum = velocity[order], spectrum[order]

        scale, offset = item["scale"], item["offset"]
        ax.plot(velocity, spectrum * scale + offset, color="black", lw=0.9, zorder=2)
        ax.axhline(offset, color="0.84", lw=0.55, zorder=0)

        broad_components = read_broad_fit_components(fitfile, sigma_min)
        for row in broad_components:
            component = broad_gaussian(velocity, row["amp"], row["v0"], row["sigma"])
            # Do not draw a colored dashed baseline across the full panel.
            visible = np.abs(component) >= 0.005 * abs(row["amp"])
            y_component = np.full_like(component, np.nan)
            y_component[visible] = component[visible] * scale + offset
            ax.plot(velocity, y_component, ls="--", lw=1.55, color=item["color"], zorder=3)
            selected_rows.append({
                "line": key, "component": row["component"],
                "v0_km_s": row["v0"], "sigma_km_s": row["sigma"],
                "amplitude_K": row["amp"],
            })

        label_x = vmin + 0.045 * (vmax - vmin)
        ax.text(label_x, offset + 0.55, item["label"], color="black", fontsize=13, va="bottom")

    ax.text(0.78, 0.92, source.upper(), transform=ax.transAxes, fontsize=19, fontweight="bold")
    ax.set_xlim(vmin, vmax)
    offsets = [item["offset"] for item in settings]
    ax.set_ylim(min(offsets) - 2.0, max(offsets) + 6.5)
    ax.set_xlabel(r"$V_{\rm LSR}\ ({\rm km/s})$", fontsize=16)
    ax.set_ylabel(r"$T_{\rm mb} (K)$", fontsize=16)
    ax.tick_params(direction="in", top=True, right=True, labelsize=12)

    fig.tight_layout()
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"[Saved] {outfile}")

    if selected_rows:
        print("line       comp    v0 [km/s]    sigma [km/s]    amp [K]")
        for row in selected_rows:
            print(f"{row['line']:<10s} {row['component']:>3d} {row['v0_km_s']:>12.3f} {row['sigma_km_s']:>15.3f} {row['amplitude_K']:>10.4f}")
    return fig, ax, selected_rows


fig_broad, ax_broad, broad_components = plot_stacked_broad_components()

