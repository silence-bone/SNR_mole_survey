"""从手动参考参数绘制一个源，复用 build_broad_line_products 的读谱与绘图函数。"""

import argparse
import copy
import json
import math
from pathlib import Path

import build_broad_line_products as profiles

# 可直接修改默认源和路径，也可通过命令行参数覆盖。
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = 'G16.0-0.5E1'
DEFAULT_REFERENCE = ROOT / 'output/broad_lines/profile_layout_reference.json'
DEFAULT_OUTPUT = ROOT / 'output/broad_lines/single_source'


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', nargs='?', default=DEFAULT_SOURCE,
                        help='目标名，例如 HC40NW、W44W 或 3C391NE')
    parser.add_argument('--reference', type=Path, default=DEFAULT_REFERENCE,
                        help='手动调整后的布局 JSON')
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT,
                        help='PNG、PDF 和实际参数记录的输出目录')
    return parser.parse_args()


def validate_number(value, label, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f'{label} 必须是数值。')
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f'{label} 必须是有限数值' + ('且大于 0。' if positive else '。'))


def validate_limits(limits, label, allow_auto=False):
    if not isinstance(limits, list) or len(limits) != 2:
        hint = '；自动范围请写 null 或 [null, null]，不要写 [null]' if allow_auto else ''
        raise ValueError(f'{label} 必须是 [最小值, 最大值]，实际收到 {limits!r}{hint}。')
    for value in limits:
        if allow_auto and value is None:
            continue
        validate_number(value, label)
    if all(value is not None for value in limits) and limits[0] >= limits[1]:
        raise ValueError(f'{label} 的最小值必须小于最大值。')


def read_reference(path, source):
    document = json.loads(path.read_text())
    layouts = document['target_layouts']
    if source not in layouts:
        raise ValueError(f'参考文件中没有 {source}；请先添加该目标的布局参数。')
    layout = copy.deepcopy(layouts[source])
    layout['transition_labels'] = profiles.reference_transition_labels(document)
    validate_limits(layout['xlim'], 'xlim (km/s)')
    if layout.get('ylim') is not None:
        validate_limits(layout['ylim'], 'ylim (K)', allow_auto=True)
    for family_name, offset in layout['baselines_K'].items():
        validate_number(offset, f'baselines_K.{family_name}')
    for transition, multiplier in layout['multipliers'].items():
        validate_number(multiplier, f'multipliers.{transition}', positive=True)
    return layout


def prepare_layout(layout, results):
    """保留手调参数，为新增跃迁补齐倍率与基线（K）。"""
    envelopes = profiles.raw_envelopes(results, *layout['xlim'])
    added = []
    for result in results:
        family_name = profiles.plot_family(result)
        if result['line'] not in layout['multipliers']:
            siblings = [layout['multipliers'][entry['line']] for entry in results
                        if profiles.plot_family(entry) == family_name
                        and entry['line'] in layout['multipliers']]
            layout['multipliers'][result['line']] = siblings[0] if siblings else 1
            added.append(result['line'])
    for result in results:
        family_name = profiles.plot_family(result)
        if family_name not in layout['baselines_K']:
            bottom = min((envelopes[entry['line']][0]
                          * layout['multipliers'][entry['line']]
                          + layout['baselines_K'][profiles.plot_family(entry)]
                          for entry in results
                          if profiles.plot_family(entry) in layout['baselines_K']), default=0)
            top = max(envelopes[entry['line']][1] * layout['multipliers'][entry['line']]
                      for entry in results if profiles.plot_family(entry) == family_name)
            layout['baselines_K'][family_name] = bottom - top - profiles.BASELINE_STEP_K
            if result['line'] not in added:
                added.append(result['line'])
    layout['automatically_configured_transitions'] = added

    requested_limits = layout.get('ylim')
    if requested_limits is None:
        requested_limits = [None, None]
    if added or any(value is None for value in requested_limits):
        envelopes = profiles.raw_envelopes(results, *layout['xlim'])
        lower_bounds, upper_bounds = [], []
        for result in results:
            offset = layout['baselines_K'][profiles.plot_family(result)]
            multiplier = layout['multipliers'][result['line']]
            low, high = envelopes[result['line']]
            lower_bounds.append(low * multiplier + offset)
            upper_bounds.append(high * multiplier + offset)
        lower, upper = min(lower_bounds), max(upper_bounds)
        margin = max(0.15, 0.035 * (upper - lower))
        automatic_limits = [lower - margin, upper + margin]
        layout['ylim'] = [automatic if requested is None else requested
                          for requested, automatic in zip(requested_limits, automatic_limits)]
        if added:
            layout['ylim'] = [min(layout['ylim'][0], automatic_limits[0]),
                              max(layout['ylim'][1], automatic_limits[1])]
    validate_limits(layout['ylim'], '最终 ylim (K)')
    return layout


def save_record(directory, source, reference, layout, versions, warnings, skipped):
    record = {
        'source': source,
        'reference_file': str(reference.resolve()),
        'display_formula': 'T_display = T_mb * multiplier + baseline_K',
        'layout': layout,
        'selected_fits': [entry for entry in versions if entry['source'] == source
                          and entry['line'] in layout['multipliers']],
        'unconfigured_transitions': skipped,
        'calibration_notes': warnings,
    }
    path = directory / f'{source}_layout_used.json'
    path.write_text(json.dumps(record, ensure_ascii=False, indent=2) + '\n')


def main():
    arguments = parse_arguments()
    source = arguments.source
    if source.upper() in profiles.W44_REGIONS:
        source = source.upper()
    layout = read_reference(arguments.reference, source)
    selected, versions = profiles.select_fit_files()
    selected = [entry for entry in selected if entry[0] == source]
    results, _, warnings = profiles.collect_results(
        selected, sources=[source], transitions=set(layout['multipliers']))
    if not results:
        raise ValueError(f'没有找到 {source} 的原始谱数据或拟合文件。')
    missing = set(layout['multipliers']) - {entry['line'] for entry in results}
    if missing:
        raise ValueError(f'参考文件中的跃迁缺少原始谱数据：{sorted(missing)}')
    for result in results:
        if result['p'] is None:
            warnings.append(f"{result['line']}: 仅绘制原始谱；input_file={result['meta']['input_file']}; "
                            f"eta={result['meta']['applied_beam_efficiency']}; 未扣基线、未合并通道。")
    layout = prepare_layout(layout, results)
    profiles.write_profiles(results, [source], {source: layout}, arguments.output_dir)
    save_record(arguments.output_dir, source, arguments.reference, layout, versions, warnings, [])
    # print(f'{source}: {arguments.output_dir / (source + ".pdf")}')
    print(f'{source}: {arguments.output_dir / (source + ".png")}')
    print('实际使用的参数已保存；参考文件保持原样。')


if __name__ == '__main__':
    main()
