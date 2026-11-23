"""Append HCO+ anchored integrated-intensity ratios to the broad-line document."""
from pathlib import Path
import re
import math
import json

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output/broad_lines'
LINES = ['12CO10', '12CO21', '12CO43', 'HCO+32', 'HCN10', 'HCN32', 'SiO21']
PAIRS = [('HCO+10', '12CO10'), ('HCO+10', '12CO21'),
         ('12CO21', '12CO10'), ('12CO43', '12CO10'), ('12CO43', '12CO21'),
         ('HCO+10', 'HCN10'), ('HCO+10', '12CO43')]
BEGIN = '% BEGIN HCO ANCHORED RATIOS'
END = '% END HCO ANCHORED RATIOS'


def main():
    versions = {}
    for path in (ROOT / 'output/line_fit').glob('*fit*.txt'):
        match = re.fullmatch(r'(.+)_(12CO\d+|HCN\d+|HCO[+p]\d+|SiO\d+)_fit(\d*)\.txt', path.name)
        if match:
            source, line, version = match.groups()
            source = source.upper() if source.lower() in {'w44w', 'w44n', 'w44c', 'w44s'} else source
            versions.setdefault((source, line.replace('HCOp', 'HCO+')), []).append((int(version or 1), path.stat().st_mtime_ns, path))
    fits = {}
    for key, choices in versions.items():
        path = max(choices)[2]
        text = path.read_text()
        legacy = 'amp_K' not in text
        components = []
        for row in text.splitlines():
            if not row.strip() or row.startswith('#'):
                continue
            a = list(map(float, row.split()))
            if legacy:
                v, ev, sigma, esigma, amp, eamp, integ, ei = a[2:10]
                width = math.sqrt(8 * math.log(2)) * sigma
            else:
                amp, eamp, v, ev, sigma, esigma, width, ew, integ, ei = a[2:12]
            if amp > 0 and a[1] > 0:
                components.append(dict(component=int(a[0]), v=v, fwhm=width, integrated=integ,
                                       error=ei, fit_file=str(path.relative_to(ROOT)), legacy=legacy))
        fits[key] = components
    audit = []
    for (source, line), components in sorted(fits.items()):
        if line != 'HCO+10':
            continue
        for anchor in sorted(components, key=lambda c: c['v']):
            if anchor['fwhm'] <= 7:
                continue
            matches = {'HCO+10': anchor}
            for other in LINES:
                candidates = [c for c in fits.get((source, other), [])
                              if c['fwhm'] > 8 and abs(c['v'] - anchor['v']) <= 10]
                matches[other] = min(candidates, key=lambda c: (abs(c['v'] - anchor['v']), c['v'], c['component'])) if candidates else None
            ratios = []
            for numerator, denominator in PAIRS:
                n, d = matches[numerator], matches[denominator]
                ratios.append(n['integrated'] / d['integrated'] if n and d and d['integrated'] > 0 else None)
            audit.append(dict(source=source, anchor=anchor, matches=matches, ratios=ratios))
    head = r'Target & $v_0$ & FWHM & $H_{10}/C_{10}$ & $H_{10}/C_{21}$ & $C_{21}/C_{10}$ & $C_{43}/C_{10}$ & $C_{43}/C_{21}$ & $H_{10}/N_{10}$ & $H_{10}/C_{43}$ \\'
    tex = [BEGIN, r'\clearpage', r'\begingroup\footnotesize', r'\setlength{\tabcolsep}{3pt}',
           r'\begin{longtable}{lrrrrrrrrr}',
           r'\caption{Integrated-intensity ratios of broad components matched to HCO$^+$ (1--0).}\label{tab:broad-ratios}\\',
           r'\toprule', head, r' & (km s$^{-1}$) & (km s$^{-1}$) & & & & & & & \\',
           r'\midrule\endfirsthead', r'\toprule', head, r'\midrule\endhead', r'\bottomrule\endfoot']
    for row in audit:
        a = row['anchor']
        tex.append(' & '.join([row['source'].replace('_', r'\_'), f"{a['v']:.2f}", f"{a['fwhm']:.2f}"] +
                              ['-' if r is None else f'{r:.3f}' for r in row['ratios']]) + r' \\')
    tex += [r'\end{longtable}',
            r'\noindent\footnotesize Notes. Each row is a positive-emission HCO$^+$ (1--0) component with FWHM $>7$ km s$^{-1}$. Counterparts in $^{12}$CO (1--0), (2--1), (4--3), HCO$^+$ (3--2), HCN (1--0), (3--2), and SiO (2--1) require FWHM $>8$ km s$^{-1}$ and $|v_0-v_{0,\mathrm{HCO^+}}|\leq10$ km s$^{-1}$. If several qualify, the nearest velocity is selected (ties: lower velocity, then component number). Matching is performed independently for each reference. $H_{10}$, $C_{J,J-1}$, and $N_{10}$ denote the fitted component integrals $\int T\,dv$ in K km s$^{-1}$ for HCO$^+$, $^{12}$CO, and HCN (1--0), respectively. A dash indicates an unavailable qualifying component. Ratios use the existing fitted temperature scales without additional beam or calibration corrections; the calibration caveats of the preceding table also apply. The highest numbered fit version is used, with modification time breaking ties.',
            r'\endgroup', END]
    fragment = '\n'.join(tex) + '\n'
    (OUT / 'broad_line_ratios.tex').write_text(fragment)
    (OUT / 'broad_line_ratios.json').write_text(json.dumps(audit, indent=2))
    document = OUT / 'broad_line_table_standalone.tex'
    original = document.read_text()
    if BEGIN in original:
        start = original.index(BEGIN)
        stop = original.index(END, start) + len(END)
        updated = original[:start] + fragment.rstrip() + original[stop:]
    else:
        updated = original.replace(r'\end{document}', fragment + '\n' + r'\end{document}')
    document.write_text(updated)
    print(f'{len(audit)} components in {len({r["source"] for r in audit})} sources')
    print('Multiple qualifying counterparts:', sum(sum(c['fwhm'] > 8 and abs(c['v']-r['anchor']['v']) <= 10 for c in fits.get((r['source'], line), [])) > 1 for r in audit for line in LINES))


if __name__ == '__main__':
    main()
