"""Rebuild IRAM-dependent plots from corrected fits, retaining W44 and MCMC runs."""
from pathlib import Path
import ast
import csv
import json
import os
import re
import shutil
import sys

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output/iram_beff_correction'
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'tmp/iram_beff_correction/mplconfig'))
sys.path.insert(0, str(ROOT / 'code/other'))
sys.path.insert(0, str(ROOT / 'code/fitting'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import build_broad_line_products as profiles
import plot_single_line_profile as single
import plot_single_table_figure as tables
import update_grouped_broad_table as groupedTable
from line_fit_core import fit_baseline
from recalibrate_iram_products import VERSION, metadata, correctedText, resolveInput
from beam_efficiency import iram_header_efficiency

generated = []


def backup(path):
    destination = OUT / 'before' / path.relative_to(ROOT)
    if path.exists() and not destination.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    generated.append(str(path.relative_to(ROOT)))


def saveFigure(figure, path):
    backup(path)
    figure.savefig(path, dpi=180, bbox_inches='tight')


def fitData(path):
    text = path.read_text()
    meta = metadata(text)
    staged = '# component role ' in text
    if not staged:
        meta, _, components = profiles.read_fit_file(path)
    else:
        components = []
        for line in text.splitlines():
            if not line.strip() or line.startswith('#'):
                continue
            row = line.split()
            components.append(dict(component=int(row[0]), amp=float(row[4]),
                                   v=float(row[6]), sigma=float(row[8]), fwhm=float(row[10])))
    velocity, temperature = profiles.load_spectrum(
        resolveInput(meta['input_file']), float(meta['applied_beam_efficiency']), meta)
    if staged and 'baseline_coefficients_in_increasing_power' not in meta:
        baseline, _, _, _ = fit_baseline(velocity, temperature,
            ast.literal_eval(meta['baseline_windows_km_s']), int(meta['baseline_order']), 3)
        temperature -= baseline
    return meta, components, velocity, temperature


def drawFits(rows):
    for row in rows:
        path = ROOT / row['path']
        meta, components, velocity, temperature = fitData(path)
        model = sum(c['amp'] * np.exp(-0.5*((velocity-c['v'])/c['sigma'])**2)
                    for c in components)
        if 'velocity_window_km_s' in meta:
            limits = ast.literal_eval(meta['velocity_window_km_s'])
        else:
            limits = (min(c['v']-3*c['sigma'] for c in components)-5,
                      max(c['v']+3*c['sigma'] for c in components)+5)
        mask = (velocity >= limits[0]) & (velocity <= limits[1])
        fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                                 gridspec_kw={'height_ratios':[3,1]}, layout='constrained')
        axes[0].step(velocity[mask], temperature[mask], where='mid', color='black', lw=.7, label='Spectrum')
        axes[0].plot(velocity[mask], model[mask], color='crimson', lw=1.3, label='Total fit')
        for c in components:
            curve = c['amp']*np.exp(-.5*((velocity[mask]-c['v'])/c['sigma'])**2)
            axes[0].plot(velocity[mask], curve, '--', lw=.7, alpha=.7)
        axes[0].set(title=path.stem, ylabel=r'$T_{mb}$ (K)', xlim=limits)
        axes[0].legend()
        axes[1].step(velocity[mask], (temperature-model)[mask], where='mid', color='.35', lw=.7)
        rms = float(meta['baseline_rms_K'])
        axes[1].axhspan(-rms, rms, color='.5', alpha=.2)
        axes[1].set(xlabel=r'$v_{LSR}$ (km s$^{-1}$)', ylabel='Residual (K)')
        target = path.with_suffix('.png')
        if not target.exists() and path.parts[len(ROOT.parts)] == 'code':
            target = OUT / 'fit_figures' / target.name
        saveFigure(fig, target)
        plt.close(fig)


def drawProfiles(results, sources):
    summaryPath = ROOT / 'output/broad_lines/summary.json'
    summary = json.loads(summaryPath.read_text())
    reference = ROOT / 'output/broad_lines/profile_layout_reference.json'
    for source in sources:
        data = [r for r in results if r['source'] == source]
        layout = single.prepare_layout(single.read_reference(reference, source), data)
        # Preserve the existing manually selected transitions and display offsets.
        fig, ax = plt.subplots(figsize=profiles.FIGURE_SIZE)
        for entry in data:
            if entry['line'] in layout['multipliers']:
                profiles.plot_transition(ax, entry, layout)
        profiles.decorate_profile(fig, ax, source, layout)
        for folder in ('figures', 'single_source'):
            for extension in ('.png', '.pdf'):
                target = ROOT / 'output/broad_lines' / folder / (source + extension)
                if folder == 'figures' or target.exists():
                    saveFigure(fig, target)
        plt.close(fig)
        summary['target_layouts'][source] = layout
        tableRows, _ = tables.load_source_components(source, reference)
        _, grouped = tables.render_table([r for r in tableRows if r['line'] in groupedTable.ORDER])
        if grouped:
            for folder in ('output/broad_lines/component_tables', 'output/broad_lines/single_source',
                           'output/fit_diagnosis_20261014/component_tables'):
                directory = ROOT / folder
                if not (directory / f'{source}_components.png').exists():
                    continue
                for ext in ('.png','.pdf'):
                    backup(directory / f'{source}_components{ext}')
                tables.draw_table(source, grouped, directory, 220)
    backup(summaryPath)
    summaryPath.write_text(json.dumps(summary, indent=2)+'\n')


def scaleWindow(record):
    """Keep each window's original input, baseline and selection; scale IRAM only."""
    for line, entry in record['lines'].items():
        meta = entry['input']
        if meta.get('calibration_revision') == VERSION:
            continue
        path = resolveInput(meta['input_file'])
        h = fits.getheader(path)
        factor = iram_header_efficiency(h, path)
        if factor == 1:
            continue
        for key in ('integrated_K_kms','error_K_kms','rms_K'):
            if entry.get(key) is not None:
                entry[key] *= factor
        text = '\n'.join(f'# {key}: {value}' for key,value in meta.items())
        updated, _ = correctedText(text, h)
        entry['input'] = metadata(updated)
    if 'ratios' in record:
        denominator = record['lines'].get('12CO43', {}).get('integrated_K_kms', 0)
        record['ratios'] = {line: entry['integrated_K_kms']/denominator
                            for line,entry in record['lines'].items() if denominator}


def drawWindow(record, directory):
    source = record['source']
    reference = single.read_reference(ROOT/'output/broad_lines/profile_layout_reference.json', source)
    fig, ax = plt.subplots(figsize=(10, 7))
    for line, entry in record['lines'].items():
        meta = entry['input']
        velocity, intensity = profiles.load_spectrum(resolveInput(meta['input_file']),
                                      float(meta.get('applied_beam_efficiency', 1)), meta)
        scale, offset = entry.get('scale',1), entry.get('offset_K',0)
        show = (velocity >= reference['xlim'][0]) & (velocity <= reference['xlim'][1])
        y = intensity*scale+offset
        ax.step(velocity[show], y[show], where='mid', lw=.8, color=profiles.COLOR.get(line),
                label=profiles.LABEL.get(line,line)+f' ×{scale:g}')
        ax.axhline(offset, color='.8', lw=.5)
        use = show & (velocity >= record['window_kms'][0]) & (velocity <= record['window_kms'][1])
        ax.fill_between(velocity[use], offset, y[use], step='mid', alpha=.25, color=profiles.COLOR.get(line))
    ax.set(title=source, xlim=reference['xlim'], xlabel=r'$v_{LSR}$ (km s$^{-1}$)',
           ylabel=r'$T_{mb}$ × scale + offset (K)')
    ax.legend(fontsize=8)
    for ext in ('.png','.pdf'):
        saveFigure(fig, directory/f'{source}_window_profiles{ext}')
    plt.close(fig)


def updateWindows():
    for path in sorted((ROOT/'output/window_profiles').glob('*_window.json')):
        if 'W44' in path.name.upper():
            continue
        record = json.loads(path.read_text())
        scaleWindow(record)
        backup(path)
        path.write_text(json.dumps(record, indent=2)+'\n')
        drawWindow(record, path.parent)


def drawDiagnostics(results, sources):
    audit = []
    for source in sources:
        path = ROOT/'output/fit_diagnosis_20261014'/f'{source}.png'
        if not path.exists():
            continue
        data = [r for r in results if r['source'] == source]
        anchors = [c for r in data if r['line']=='HCO+10' for c in r['comps'] if c['amp']>0 and c['fwhm']>7]
        low = max(-40, min(c['v']-c['fwhm'] for c in anchors)-10)
        high = min(160, max(c['v']+c['fwhm'] for c in anchors)+10)
        fig, axes = plt.subplots(4,1,figsize=(11,10),sharex=True,layout='constrained')
        for ax,line in zip(axes, ['HCO+10','12CO10','12CO21','12CO43']):
            entry = next((r for r in data if r['line']==line), None)
            if entry is None:
                continue
            x,y = entry['x'],entry['y']
            use = (x>=low)&(x<=high)
            model = np.zeros_like(x)
            for c in entry['comps']:
                curve = c['amp']*np.exp(-.5*((x-c['v'])/c['sigma'])**2)
                model += curve
                ax.plot(x[use],curve[use],'--',lw=.6)
                audit.append(dict(source=source,line=line,fit_file=entry['p'].name,**c,
                                  rms=entry['meta'].get('baseline_rms_K'),
                                  input_file=entry['meta'].get('input_file')))
            ax.step(x[use],y[use],where='mid',color='k',lw=.6)
            ax.plot(x[use],model[use],color='crimson',lw=1)
            ax.plot(x[use],(y-model)[use],color='.5',lw=.5)
            ax.set(ylabel=line+'\nT (K)',xlim=(low,high))
        axes[0].set_title(source)
        axes[-1].set_xlabel('LSR velocity (km/s)')
        saveFigure(fig,path)
        plt.close(fig)
    path=ROOT/'output/fit_diagnosis_20261014/components.csv'
    backup(path)
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(audit[0]));writer.writeheader();writer.writerows(audit)


def main():
    rows = list(csv.DictReader((OUT/'corrected_fits.csv').open()))
    if '--skip-fits' not in sys.argv:
        drawFits(rows)
    selected, versions = profiles.select_fit_files()
    results, broadRows, warnings = profiles.collect_results(
        [entry for entry in selected if 'W44' not in entry[0].upper()])
    oldRows = list(csv.DictReader((ROOT/'output/broad_lines/broad_components.csv').open()))
    broadRows.extend(row for row in oldRows if 'W44' in row['source'].upper())
    oldManifest = json.loads((ROOT/'output/broad_lines/selection_manifest.json').read_text())
    warnings.extend(w for w in oldManifest.get('warnings',[]) if 'W44' in w.upper())
    changed = {Path(row['path']).name for row in rows if row['path'].startswith('code/')}
    sources = sorted({r['source'] for r in results if r['p'] and r['p'].name in changed})
    for name in ('broad_components.csv','selection_manifest.json'):
        backup(ROOT/'output/broad_lines'/name)
    profiles.write_inventory(broadRows,versions,warnings)
    drawProfiles(results,sources)
    updateWindows()
    drawDiagnostics(results,sources)
    (OUT/'regenerated_files.json').write_text(json.dumps(sorted(set(generated)),indent=2)+'\n')
    print(f'Regenerated {len(set(generated))} outputs for {len(sources)} non-W44 targets.')


if __name__ == '__main__':
    main()
