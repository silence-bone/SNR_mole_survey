"""Plot offset CO/HCO+/HCN spectra and measure unscaled window integrals."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import build_broad_line_products as profiles
from plot_single_line_profile import read_reference, DEFAULT_REFERENCE
from profile_plot_style import FIGURE_SIZE, label_spectrum, plot_context, style_axes

ROOT = Path(__file__).resolve().parents[2]
# User settings: edit these defaults here, or override target/window via CLI.
TARGET = 'G16.0-0.5E2'
RATIO_WINDOW = (70.0, 80.0)  # km/s
REFERENCE_JSON = DEFAULT_REFERENCE
OUTPUT_DIR = ROOT / 'code/window_profiles'
LINES = ['12CO10', '12CO21', '12CO43', 'HCO+10', 'HCO+32', 'HCN10', 'HCN32']
OFFSETS_K = dict(zip(LINES, [15, 10, 5, 0, -5, -10, -15, -20]))
RMS_OVERRIDES_K = {}  # Optional {line: RMS in K}, e.g. for spectra without fit RMS.
FIGSIZE = FIGURE_SIZE
DPI = 200
SHADE_ALPHA = 0.3
YLIM_PADDING_K = (1, 2)


def integrate_window(x, y, low, high):
    """Integrate channel overlap, including fractional boundary channels."""
    order = np.argsort(x)
    x, y = np.asarray(x)[order], np.asarray(y)[order]
    edges = np.r_[x[0] - (x[1]-x[0])/2, (x[:-1]+x[1:])/2,
                  x[-1] + (x[-1]-x[-2])/2]
    weights = np.maximum(0, np.minimum(edges[1:], high)-np.maximum(edges[:-1], low))
    use = weights > 0
    if low < edges[0] or high > edges[-1] or not np.all(np.isfinite(y[use])):
        raise ValueError('Window is not fully covered by finite channels')
    return float(weights[use] @ y[use]), float(np.sqrt(np.sum(weights**2)))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', nargs='?', default=TARGET)
    parser.add_argument('--ratio-window', '--window', dest='window', nargs=2,
                        type=float, default=RATIO_WINDOW)
    parser.add_argument('--reference', type=Path, default=REFERENCE_JSON)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    if args.window[0] >= args.window[1]:
        parser.error('window must be increasing')
    return args


def load_target(args):
    """Load the selected target and its reference scaling without target-specific fits."""
    source = args.source.upper() if args.source.upper() in profiles.W44_REGIONS else args.source
    layout = read_reference(args.reference, source)
    layout['xlim'] = [min(layout['xlim'][0], args.window[0]),
                      max(layout['xlim'][1], args.window[1])]
    selected, _ = profiles.select_fit_files()
    results, _, warnings = profiles.collect_results(
        [row for row in selected if row[0] == source], sources=[source], transitions=set(LINES))
    if not results:
        raise ValueError(f'No CO/HCO+/HCN spectra found for {source}')
    for line in LINES:
        layout['multipliers'].setdefault(line, 1)
    return source, layout, results, warnings


def measure_lines(source, results, layout, window, warnings):
    """Measure unscaled spectra; missing RMS is recorded without inventing errors."""
    record = dict(source=source, window_kms=list(window), warnings=warnings, lines={})
    for result in sorted(results, key=lambda row: LINES.index(row['line'])):
        line = result['line']
        scale = layout['multipliers'][line]
        value, noise_weight = integrate_window(result['x'], result['y'], *window)
        rms_value = RMS_OVERRIDES_K.get(line, result['meta'].get('baseline_rms_K'))
        rms = float(rms_value) if rms_value is not None else None
        if rms is not None and (not np.isfinite(rms) or rms <= 0):
            raise ValueError(f'Invalid RMS for {line}')
        if rms is None:
            warnings.append(f'{line}: no RMS; integral error unavailable, set RMS_OVERRIDES_K before MCMC')
        record['lines'][line] = dict(
            integrated_K_kms=value, error_K_kms=None if rms is None else rms*noise_weight,
            rms_K=rms, rms_source='override' if line in RMS_OVERRIDES_K else str(result['p']),
            offset_K=OFFSETS_K[line], rms_to_error_kms=noise_weight,
            scale=scale, input=result['meta'])
    values = record['lines']
    denominator = values.get('12CO43', {}).get('integrated_K_kms', 0)
    record['ratios'] = {line: values[line]['integrated_K_kms']/denominator
                        for line in LINES if line in values and denominator != 0}
    return record


@plot_context()
def plot_profiles(source, results, layout, window, output_dir):
    """Fill only the window's channel areas between each spectrum and its baseline."""
    fig, ax = plt.subplots(figsize=FIGSIZE)
    offsets = OFFSETS_K
    for result in sorted(results, key=lambda row: LINES.index(row['line'])):
        line = result['line']
        scale = layout['multipliers'][line]
        x, y = result['x'], result['y']*scale + offsets[line]
        ax.step(x, y, where='mid', color=profiles.COLOR[line])
        ax.axhline(offsets[line], color='0.7', lw=0.6)
        label_spectrum(ax, layout.get('transition_labels', profiles.LABEL)[line],
                       scale, offsets[line])
        edges = np.r_[x[0]-(x[1]-x[0])/2, (x[:-1]+x[1:])/2,
                      x[-1]+(x[-1]-x[-2])/2]
        left = np.maximum(edges[:-1], window[0])
        right = np.minimum(edges[1:], window[1])
        use = right > left
        xx = np.column_stack([left[use], right[use]]).ravel()
        yy = np.repeat(y[use], 2)
        ax.fill_between(xx, offsets[line], yy, color=profiles.COLOR[line], alpha=SHADE_ALPHA)
    visible = np.concatenate([
        (r['y']*layout['multipliers'][r['line']]+offsets[r['line']])[
            (r['x'] >= layout['xlim'][0]) & (r['x'] <= layout['xlim'][1])]
        for r in results])
    ax.set_ylim(np.nanmin(visible)-YLIM_PADDING_K[0], np.nanmax(visible)+YLIM_PADDING_K[1])
    ax.set(xlim=layout['xlim'], xlabel=r'$v_{LSR}$ (km s$^{-1}$)',
           ylabel=r'$T_{mb}$ × scale + offset (K)', title=source)
    style_axes(ax, source)
    fig.tight_layout()
    fig.savefig(output_dir/f'{source}_window_profiles.png', dpi=DPI)
    plt.close(fig)


def main():
    args = parse_args()
    source, layout, results, warnings = load_target(args)
    record = measure_lines(source, results, layout, args.window, warnings)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_profiles(source, results, layout, args.window, args.output_dir)
    (args.output_dir/f'{source}_window.json').write_text(json.dumps(record, indent=2))
    print("save figure to ", args.output_dir/f'{source}_window_profiles.png')
    for warning in warnings:
        print(warning)


if __name__ == '__main__':
    main()
