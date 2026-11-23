"""绘制单个源的宽线分量表图，与单源 line profile 使用相同的跃迁参考。"""

import argparse
import csv
import json
from pathlib import Path

import build_broad_line_products as profiles
from update_grouped_broad_table import render_table

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = 'HC40NW'
DEFAULT_REFERENCE = ROOT / 'output/broad_lines/profile_layout_reference.json'
DEFAULT_OUTPUT = ROOT / 'output/broad_lines/single_source'
PNG_DPI = 220
FIGURE_WIDTH_INCH = 12
FONT_SIZE = 12
FONT_FAMILY = 'Times New Roman'
ROW_HEIGHT_INCH = 0.28
HEADERS = ['Target', 'Transition', r'$v_0$', r'$\mathrm{T}_p$', 'FWHM', r'$\int T\,dv$']
UNITS = ['', '', r'km s$^{-1}$', 'K', r'km s$^{-1}$', r'K km s$^{-1}$']
COLUMN_WIDTHS = [0.15, 0.12, 0.13, 0.13, 0.13, 0.13]
PARAMETERS = [('v', 'ev'), ('amp', 'eamp'), ('fwhm', 'efwhm'),
              ('integrated', 'eintegrated')]


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', nargs='?', default=DEFAULT_SOURCE)
    parser.add_argument('--reference', type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--dpi', type=int, default=PNG_DPI,
                        help='直接渲染 PNG 的分辨率；PDF 始终为矢量输出')
    arguments = parser.parse_args()
    if arguments.dpi <= 0:
        parser.error('--dpi 必须大于 0')
    return arguments


def load_source_components(source, reference):
    """仅需拟合参数，无需重新加载 FITS 或对强度应用显示倍率。"""
    selected, versions = profiles.select_fit_files()
    selected = [entry for entry in selected if entry[0] == source]
    if not selected:
        raise ValueError(f'没有找到 {source} 的拟合文件。')
    rows = []
    for _, transition, path in selected:
        _, legacy, components = profiles.read_fit_file(path)
        for component in components:
            if component['type'] > 0:
                rows.append(dict(source=source, line=transition, fit_file=path.name,
                                 legacy=legacy, **component))
    versions = [entry for entry in versions if entry['source'] == source]
    return rows, versions


def table_cells(source, grouped, labels=None):
    """沿用已有 HCO+ 速度分组，表中强度均为拟合原值。"""
    labels = profiles.LABEL if labels is None else labels
    cells = [HEADERS, UNITS]
    section_rows = []
    previous_group = None
    has_reference = any(row['matched'] for row in grouped)
    for row in grouped:
        first_in_group = row['group'] != previous_group
        if first_in_group:
            if row['matched']:
                title = ('Reference ' + labels['HCO+10'] + r': $v_0=$ '
                         + f'{float(row["reference_v0"]):.2f}' + r' km s$^{-1}$')
            else:
                title = ('Unmatched components' if has_reference else
                         'No broad ' + labels['HCO+10'] + ' reference')
            section_rows.append((len(cells), title))
            cells.append([''] * len(HEADERS))
        cells.append(['', labels[row['line']]] + [
            profiles.fmt(float(row[value]), float(row[error])) for value, error in PARAMETERS
        ])
        previous_group = row['group']
    return cells, section_rows


def draw_table(source, grouped, directory, dpi, labels=None):
    cells, section_rows = table_cells(source, grouped, labels)
    row_count = len(cells)
    table_height = ROW_HEIGHT_INCH * row_count
    profiles.plt.rcParams.update({
        'font.family': FONT_FAMILY, 'font.size': FONT_SIZE, 'pdf.fonttype': 42,
        'mathtext.fontset': 'custom',
        'mathtext.rm': FONT_FAMILY,
        'mathtext.it': f'{FONT_FAMILY}:italic',
        'mathtext.bf': f'{FONT_FAMILY}:bold',
        'mathtext.cal': f'{FONT_FAMILY}:italic',
    })
    figure = profiles.plt.figure(figsize=(FIGURE_WIDTH_INCH, table_height + 1.2))
    axes = figure.add_axes([0.03, 0.65 / (table_height + 1.2), 0.94,
                            table_height / (table_height + 1.2)])
    axes.axis('off')
    table = axes.table(cellText=cells, colWidths=COLUMN_WIDTHS,
                       cellLoc='center', bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(FONT_SIZE)
    for (row_index, _), cell in table.get_celld().items():
        cell.visible_edges = ''
        if row_index == 0:
            cell.visible_edges = 'T'
        if row_index in (1, row_count - 1):
            cell.visible_edges += 'B'
        cell.set_linewidth(0.7)
    for row_index, title in section_rows:
        # 分组线与标题仅覆盖后五列，Target 列保持为完整的跨行区域。
        if row_index > 2:
            for column in range(1, len(HEADERS)):
                table[row_index, column].visible_edges = 'T'
                table[row_index, column].set_linewidth(0.4)
        axes.text((1 + COLUMN_WIDTHS[0]) / 2, 1 - (row_index + 0.5) / row_count, title,
                  transform=axes.transAxes, ha='center', va='center', fontsize=FONT_SIZE)
    target = source + (r'$^{\dagger}$' if any(row['legacy'] for row in grouped) else '')
    axes.text(COLUMN_WIDTHS[0] / 2, (1 - 2 / row_count) / 2, target,
              transform=axes.transAxes, ha='center', va='center', fontsize=FONT_SIZE)
    axes.set_title(f'{source}: broad emission components', fontsize=12, pad=12)
    # note = 'Positive emission: HCO+ (1–0) FWHM > 7 km/s; other lines > 8 km/s. Matching: |Δv| ≤ 10 km/s.'
    # if any(row['legacy'] for row in grouped):
    #     note += '\n† Legacy calibration metadata unavailable; fitted values retained.'
    # figure.text(0.5, 0.12 / (table_height + 1.2), note, ha='center', va='bottom', fontsize=8)
    # 从文字、线条直接输出；不经过 JPEG，也不放大低分辨率位图。
    figure.savefig(directory / f'{source}_components.png',
                   dpi=dpi, bbox_inches='tight', pad_inches=0.12)
    figure.savefig(directory / f'{source}_components.pdf', bbox_inches='tight', pad_inches=0.12)
    profiles.plt.close(figure)


def save_records(source, grouped, versions, arguments):
    path = arguments.output_dir / f'{source}_components.csv'
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(grouped[0]))
        writer.writeheader()
        writer.writerows(grouped)
    record = {
        'source': source, 'reference_file': str(arguments.reference.resolve()),
        'transition_labels': profiles.reference_transition_labels(
            json.loads(arguments.reference.read_text())),
        'selected_fits': versions, 'display_rows': len(grouped), 'png_dpi': arguments.dpi,
        'selection': 'Positive emission; HCO+10 FWHM > 7 km/s, others > 8 km/s.',
        'grouping': 'Existing HCO+ reference matching within 10 km/s; unmatched retained.',
        'intensity': 'Original fitted values; profile offsets and multipliers are not applied.',
    }
    (arguments.output_dir / f'{source}_table_used.json').write_text(
        json.dumps(record, ensure_ascii=False, indent=2) + '\n'
    )


def main():
    arguments = parse_arguments()
    source = arguments.source
    if source.upper() in profiles.W44_REGIONS:
        source = source.upper()
    rows, versions = load_source_components(source, arguments.reference)
    _, grouped = render_table(rows)
    if not grouped:
        raise ValueError(f'{source} 的所选跃迁没有满足阈值的宽发射分量。')
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    labels = profiles.reference_transition_labels(json.loads(arguments.reference.read_text()))
    draw_table(source, grouped, arguments.output_dir, arguments.dpi, labels)
    # save_records(source, grouped, versions, arguments)
    print(f'{source}: {len(grouped)} rows; {arguments.output_dir / (source + "_components.png")}')



if __name__ == '__main__':
    main()
