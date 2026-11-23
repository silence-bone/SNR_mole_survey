"""汇总宽线拟合结果，并生成按分子分层显示的 line profile。

运行时保留已有拟合值与温标处理；这里只读取、整理和绘图，不重新拟合。
使用 --figures-only 可跳过 LaTeX 和 README 的生成，但仍会更新分量清单。
"""

import csv
import json
import math
import os
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'code' / 'other'))
OUT = ROOT / 'output' / 'broad_lines'
FIT_DIRECTORY = ROOT / 'code' / 'fitting' / 'line_fit'
FIGURE_DIRECTORY = OUT / 'figures'

# 用户可调整的参数：速度单位 km/s，温度与基线间距单位 K。
BROAD_FWHM_KMS = 8
SPEED_OF_LIGHT_KMS = 299792.458
MINIMUM_FAMILY_GAP_K = 0.3
BASELINE_STEP_K = 5
MINIMUM_MULTIPLIER = 10
MAXIMUM_MULTIPLIER = 50
TARGET_FILL_FRACTION = 0.75
FIGURE_SIZE = (8.0, 10.0)
PNG_DPI = 200
BROAD_COMPONENT_LINEWIDTH = 1.8
PROFILE_FONT_SIZES = {
    'base': 14,
    'title': 20,
    'axis_label': 20,
    'tick': 14,
    'legend': 13,
    'baseline': 11,
    'caption': 11,
}
LEGEND_LOCATION = 'upper left'
LEGEND_AUTO_SIDE = True
AUTO_LEGEND_HEADROOM = True
LEGEND_CLEARANCE_POINTS = 6
VELOCITY_OVERRIDES = {'3C391SW': (-20, 180)}
W44_REGIONS = {'W44W', 'W44N', 'W44C', 'W44S'}
APEX_FAMILY = 'apex CO(2–1)'
ARCHIVE_FAMILY = 'archive CO(1–0)'
FAMILIES = ['CO', 'HCO+', 'HCN', 'SiO']
ORDER = ['12CO10', '12CO21', '12CO21_apex', '12CO43', 'HCN10', 'HCN32', 'HCO+10', 'HCO+32', 'SiO21']
LEGEND_ORDER = ['12CO10', '12CO21', '12CO21_apex', '12CO43', 'HCO+10', 'HCO+32', 'HCN10', 'HCN32', 'SiO21']
LABEL = {
    '12CO10': r'$^{12}$CO (1–0)',
    '12CO21': r'$^{12}$CO (2–1)',
    '12CO21_apex': '12CO(2-1) (apex)',
    '12CO43': r'$^{12}$CO (4–3)',
    'HCN10': 'HCN (1–0)',
    'HCN32': 'HCN (3–2)',
    'HCO+10': r'HCO$^+$ (1–0)',
    'HCO+32': r'HCO$^+$ (3–2)',
    'SiO21': 'SiO (2–1)',
}
COLOR = dict(zip(ORDER, [
    '#777777', '#0072B2', '#56B4E9', '#D55E00', '#009E73',
    '#CC79A7', '#E69F00', '#6A3D9A', '#332288',
]))
# 旧版拟合文件没有完整的温标元数据，保留此前采用的主波束效率。
LEGACY_EFF = {'12CO10': 0.75, '12CO21': 0.59, '12CO43': 0.61,
              'HCO+10': 0.78, 'HCO+32': 0.82}
# 保留 summary.json 中的历史字段；实际绘图倍率由每个目标的 layout 决定。
FACTOR = {
    transition: 1 if transition.startswith('12CO') else
    30 if transition.startswith('SiO') else 10
    for transition in ORDER
}
FIT_NAME_PATTERN = re.compile(
    r'(.+)_(12CO\d+|HCN\d+|HCO[+p]\d+|SiO\d+)(_apex)?_fit(\d*)\.txt$'
)

# 在导入 matplotlib 前指定缓存位置；无界面后端适合批量生成图片。
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'tmp' / 'mplconfig'))
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from beam_efficiency import resolve_beam_efficiency, temperature_divisor, iram_header_efficiency

SIGMA_TO_FWHM = 2 * np.sqrt(2 * np.log(2))


def select_fit_files():
    """每个目标、跃迁选最高版本；同版本按修改时间选择最新文件。"""
    grouped = {}
    for path in sorted(FIT_DIRECTORY.glob('*fit*.txt')):
        match = FIT_NAME_PATTERN.match(path.name)
        if not match:
            raise ValueError(path)
        source, transition, telescope_suffix, version = match.groups()
        transition += telescope_suffix or ''
        source = source.upper() if source.upper() in W44_REGIONS else source
        transition = transition.replace('HCOp', 'HCO+')
        grouped.setdefault((source, transition), []).append((int(version or 1), path))

    selected = []
    versions = []
    for (source, transition), candidates in sorted(grouped.items()):
        _, path = max(candidates, key=lambda item: (item[0], item[1].stat().st_mtime_ns))
        selected.append((source, transition, path))
        versions.append({
            'source': source, 'line': transition, 'selected': path.name,
            'alternatives': [other.name for _, other in candidates if other != path],
        })
    return selected, versions


def read_fit_file(path):
    """将新旧两种列顺序转换为同一分量字典，保留输出字段名。"""
    text = path.read_text()
    legacy = 'amp_K' not in text
    metadata = {}
    components = []
    for line in text.splitlines():
        if line.startswith('#') and ':' in line:
            key, value = line[1:].split(':', 1)
            metadata[key.strip()] = value.strip()
        if not line.strip() or line.startswith('#'):
            continue
        values = list(map(float, line.split()))
        if legacy:
            (velocity, velocity_error, sigma, sigma_error, amplitude,
             amplitude_error, integral, integral_error) = values[2:10]
            width = SIGMA_TO_FWHM * sigma
            width_error = SIGMA_TO_FWHM * sigma_error
        else:
            (amplitude, amplitude_error, velocity, velocity_error, sigma, sigma_error,
             width, width_error, integral, integral_error) = values[2:12]
        assert np.isclose(width, SIGMA_TO_FWHM * sigma, rtol=1e-4)
        components.append(dict(
            component=int(values[0]), type=int(values[1]),
            amp=amplitude, eamp=amplitude_error, v=velocity, ev=velocity_error,
            sigma=sigma, fwhm=width, efwhm=width_error,
            integrated=integral, eintegrated=integral_error,
        ))
    return metadata, legacy, components


def spectrum_input(source, transition, path, metadata, legacy, warnings):
    if legacy and not metadata.get('input_file'):
        input_source = source.lower() if source in W44_REGIONS else source
        input_line = transition.replace('HCO+', 'HCOp')
        input_path = ROOT / f'{input_source}_{input_line}_re.fits'
        # W44 旧输入已移入 w44/；根目录仍存在旧文件时优先保持原引用。
        if not input_path.exists() and source in W44_REGIONS:
            input_path = ROOT / 'w44' / input_path.name
        efficiency = LEGACY_EFF[transition]
        warnings.append(
            f'{path.name}: legacy temperature calibration undocumented; observed '
            f'data divided by eta={efficiency}; fit values preserved.'
        )
    else:
        input_path = Path(metadata['input_file'])
        # 补齐元数据的旧列格式也按明确记录读取，列顺序仍由 legacy 决定。
        efficiency = float(metadata.get('applied_beam_efficiency', 1))
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    return input_path, efficiency


def load_spectrum(input_path, efficiency, metadata):
    header = fits.getheader(input_path)
    temperature = np.asarray(fits.getdata(input_path), float).squeeze()
    pixels = np.arange(temperature.size) + 1  # FITS 像素坐标从 1 开始。
    if 'DELTAV' in header and 'VELO-LSR' in header:
        velocity = (header['VELO-LSR'] +
                    (pixels - header['CRPIX1']) * header['DELTAV']) / 1000
    elif str(header.get('CTYPE1', '')).upper() == 'VRAD':
        velocity = header['CRVAL1'] + (pixels - header['CRPIX1']) * header['CDELT1']
        if str(header.get('CUNIT1', 'm/s')).strip().lower() != 'km/s':
            velocity = velocity / 1000
    else:
        frequency = header['CRVAL1'] + (pixels - header['CRPIX1']) * header['CDELT1']
        rest_frequency = header.get('RESTFRQ', header.get('RESTFREQ'))
        velocity = SPEED_OF_LIGHT_KMS * (1 - frequency / rest_frequency)
    temperature = temperature / efficiency
    indices = np.argsort(velocity)
    velocity, temperature = velocity[indices], temperature[indices]

    # 复现拟合时的通道合并：连续、不重叠地取均值，舍弃末尾不足一组的通道。
    channel_width = abs(np.median(np.diff(velocity)))
    target_width = float(metadata.get('fitted_channel_width_km_s', channel_width))
    block_size = max(1, int(round(target_width / channel_width)))
    if block_size > 1:
        stop = len(velocity) // block_size * block_size
        velocity = velocity[:stop].reshape(-1, block_size).mean(1)
        temperature = temperature[:stop].reshape(-1, block_size).mean(1)

    # 多项式系数按常数项到高阶项排列；没有多项式时扣除常数基线。
    if 'baseline_coefficients_in_increasing_power' in metadata:
        coefficients = list(map(
            float, metadata['baseline_coefficients_in_increasing_power'].split()
        ))
        temperature -= np.polynomial.polynomial.polyval(velocity, coefficients)
    else:
        temperature -= float(metadata.get('fitted_constant_baseline_K', 0))
    return velocity, temperature


def collect_raw_results(sources, existing, transitions=None):
    """从标准数据目录补充没有拟合的跃迁；重名时优先 snr_data。"""
    results = []
    pattern = re.compile(r'(.+)_(12CO\d+|HCN\d+|HCO[+p]\d+|SiO\d+)(_apex)?(?:_re)?\.fits$')
    for directory in (ROOT / 'snr_data', ROOT / 'w44', ROOT):
        for path in sorted(directory.glob('*.fits')):
            match = pattern.fullmatch(path.name)
            if not match:
                continue
            source, transition, telescope_suffix = match.groups()
            transition += telescope_suffix or ''
            source = source.upper() if source.upper() in W44_REGIONS else source
            transition = transition.replace('HCOp', 'HCO+')
            if transitions is not None and transition not in transitions:
                continue
            key = (source, transition)
            if source not in sources or transition not in ORDER or key in existing:
                continue
            header = fits.getheader(path)
            is_tmb = (('TMB' in str(header.get('BUNIT', '')).upper()
                       and iram_header_efficiency(header, path) == 1)
                      or str(header.get('TEMPSCAL', '')).upper() == 'TMB')
            efficiency = 1.0 if is_tmb else float(resolve_beam_efficiency(header)['beff'])
            if not is_tmb:
                efficiency = temperature_divisor(header, efficiency, path)
            metadata = {'input_file': str(path), 'applied_beam_efficiency': efficiency}
            velocity, temperature = load_spectrum(path, efficiency, metadata)
            results.append(dict(source=source, line=transition, p=None, meta=metadata,
                                legacy=False, x=velocity, y=temperature, comps=[]))
            existing.add(key)
    return results


def collect_results(selected, sources=None, transitions=None):
    results, rows, warnings = [], [], []
    for source, transition, path in selected:
        if transitions is not None and transition not in transitions:
            continue
        metadata, legacy, components = read_fit_file(path)
        for component in components:
            # 宽度清单同时保留发射和吸收；这里不按峰温正负筛选。
            if component['fwhm'] > BROAD_FWHM_KMS:
                rows.append(dict(source=source, line=transition, fit_file=path.name,
                                 legacy=legacy, **component))
        input_path, efficiency = spectrum_input(
            source, transition, path, metadata, legacy, warnings
        )
        velocity, temperature = load_spectrum(input_path, efficiency, metadata)
        results.append(dict(source=source, line=transition, p=path, meta=metadata,
                            legacy=legacy, x=velocity, y=temperature, comps=components))
    sources = set(sources if sources is not None else (entry[0] for entry in selected))
    results.extend(collect_raw_results(
        sources, {(r['source'], r['line']) for r in results}, transitions))
    rows.sort(key=lambda row: (row['source'], ORDER.index(row['line']), row['v']))
    return results, rows, warnings


def write_inventory(rows, versions, warnings):
    with (OUT / 'broad_components.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest = dict(threshold_km_s=BROAD_FWHM_KMS, versions=versions, warnings=warnings)
    (OUT / 'selection_manifest.json').write_text(json.dumps(manifest, indent=2))


def family(transition):
    if transition.startswith('12CO'):
        return 'CO'
    if transition.startswith('HCO'):
        return 'HCO+'
    if transition.startswith('HCN'):
        return 'HCN'
    return 'SiO'


def plot_family(result):
    if result['line'].endswith('_apex'):
        return APEX_FAMILY
    if result['line'] == '12CO10' and 'archive_data' in result['meta'].get('input_file', ''):
        return ARCHIVE_FAMILY
    return family(result['line'])


def model_at(result, velocity):
    return sum((
        component['amp'] * np.exp(-0.5 * (
            (velocity - component['v']) / component['sigma']
        ) ** 2)
        for component in result['comps']
    ), np.zeros_like(velocity))


def velocity_limits(data):
    components = [component for result in data for component in result['comps']]
    if not components:
        return VELOCITY_OVERRIDES.get(
            data[0]['source'], (min(float(r['x'].min()) for r in data),
                                max(float(r['x'].max()) for r in data)))
    left = min(component['v'] - 3 * component['sigma'] for component in components)
    right = max(component['v'] + 3 * component['sigma'] for component in components)
    padding = max(5, 0.08 * (right - left))
    limits = (10 * math.floor((left - padding) / 10),
              10 * math.ceil((right + padding) / 10))
    return VELOCITY_OVERRIDES.get(data[0]['source'], limits)


def raw_envelopes(data, lower, upper):
    """同时检查观测谱和总拟合，避免纵向平移后谱线重叠或超出坐标。"""
    envelopes = {}
    grid = np.linspace(lower, upper, 2000)
    for result in data:
        mask = ((result['x'] >= lower) & (result['x'] <= upper) &
                np.isfinite(result['y']))
        values = np.r_[result['y'][mask], model_at(result, grid)]
        envelopes[result['line']] = [min(0., float(values.min())),
                                     max(0., float(values.max()))]
    return envelopes


def choose_family_scale(source, family_name, entries, raw, upper_low):
    """整组跃迁共用倍率；在三个候选基线间距中选最接近目标填充率的一组。"""
    peak = max(raw[result['line']][1] for result in entries)
    minimum = max(BASELINE_STEP_K, BASELINE_STEP_K * math.ceil(
        (max(raw[result['line']][1] * MINIMUM_MULTIPLIER for result in entries)
         - upper_low + MINIMUM_FAMILY_GAP_K) / BASELINE_STEP_K
    ))
    cap = (raw.get('12CO21', [0., float('inf')])[1]
           if family_name in ('HCO+', 'HCN') else float('inf'))
    if peak > 0 and peak * MINIMUM_MULTIPLIER > cap:
        raise ValueError(f'{source} {family_name}: even x10 exceeds CO(2–1) peak')

    options = []
    for distance in (minimum, minimum + BASELINE_STEP_K, minimum + 2 * BASELINE_STEP_K):
        available = min(distance + upper_low - MINIMUM_FAMILY_GAP_K, cap)
        if peak <= 0:
            multiplier = MINIMUM_MULTIPLIER
        else:
            limit = min(MAXIMUM_MULTIPLIER, int(available / peak) // 10 * 10)
            candidates = range(MINIMUM_MULTIPLIER, limit + 1, 10)
            multiplier = min(candidates, key=lambda value: abs(
                value * peak / distance - TARGET_FILL_FRACTION
            ))
        # 在填充率接近时偏好较小间距；统一倍率保持同分子的线强比。
        score = (abs(multiplier * peak / distance - TARGET_FILL_FRACTION)
                 + 0.025 * (distance - minimum) / BASELINE_STEP_K)
        options.append((score, distance, {result['line']: multiplier for result in entries}))
    _, distance, chosen = min(options, key=lambda option: option[0])
    return distance, chosen


def layout_for(data):
    lower, upper = velocity_limits(data)
    raw = raw_envelopes(data, lower, upper)
    extents, scales, offsets = {}, {}, {}
    previous = None
    for family_name in FAMILIES:
        entries = [result for result in data if plot_family(result) == family_name]
        if not entries:
            continue
        if family_name == 'CO':
            offset = 0.
            scales.update({result['line']: 1 for result in entries})
        else:
            upper_low = extents[previous][0] if previous else 0.
            distance, chosen = choose_family_scale(
                data[0]['source'], family_name, entries, raw, upper_low
            )
            scales.update(chosen)
            offset = (offsets[previous] if previous else 0.) - distance
        extents[family_name] = [
            min(raw[result['line']][0] * scales[result['line']] for result in entries),
            max(raw[result['line']][1] * scales[result['line']] for result in entries),
        ]
        if family_name in ('HCO+', 'HCN') and '12CO21' in raw:
            assert extents[family_name][1] <= raw['12CO21'][1] + 1e-10
        assert len({scales[result['line']] for result in entries}) == 1
        offsets[family_name] = offset
        previous = family_name

    # 档案 CO(1–0) 单独放在最上方，倍率保持为 1。
    for special_family in (ARCHIVE_FAMILY, APEX_FAMILY):
        archive = [result for result in data if plot_family(result) == special_family]
        if archive:
            result = archive[0]
            low, high = raw[result['line']]
            scales[result['line']] = 1
            top = max(offsets[name] + extents[name][1] for name in offsets) if offsets else 0.
            offsets[special_family] = max(BASELINE_STEP_K, BASELINE_STEP_K * math.ceil(
                (top - low + MINIMUM_FAMILY_GAP_K) / BASELINE_STEP_K
            ))
            extents[special_family] = [low, high]
    top = max(extents[name][1] + offsets[name] for name in offsets)
    bottom = min(extents[name][0] + offsets[name] for name in offsets)
    margin = max(0.15, 0.035 * (top - bottom))
    ordered = sorted(offsets, key=offsets.get, reverse=True)
    for upper_family, lower_family in zip(ordered, ordered[1:]):
        assert (offsets[lower_family] + extents[lower_family][1] <
                offsets[upper_family] + extents[upper_family][0]), data[0]['source']
    return dict(xlim=[lower, upper], ylim=[bottom - margin, top + margin],
                baselines_K=offsets, envelopes_K=extents,
                gap_K=MINIMUM_FAMILY_GAP_K, multipliers=scales)


def reference_transition_labels(document):
    """全局显示名称不改变输入文件、跃迁键或拟合参数。"""
    overrides = document.get('global_rules', {}).get('transition_labels', {})
    if not isinstance(overrides, dict):
        raise ValueError('global_rules.transition_labels 必须是跃迁名称字典。')
    for transition, label in overrides.items():
        if transition not in LABEL or not isinstance(label, str) or not label.strip():
            raise ValueError(f'无效的全局跃迁标签：{transition} = {label!r}')
    return {**LABEL, **overrides}


def plot_transition(axes, result, layout):
    multiplier = layout['multipliers'][result['line']]
    offset = layout['baselines_K'][plot_family(result)]
    color = COLOR[result['line']]
    label = layout.get('transition_labels', LABEL)[result['line']]
    axes.step(result['x'], result['y'] * multiplier + offset, where='mid',
              color=color, lw=0.65, alpha=0.8, label=label + f' ×{multiplier}',
              gid=result['line'])
    if not result['comps']:
        return
    grid = np.linspace(*layout['xlim'], 1600)
    for component in result['comps']:
        if component['fwhm'] > BROAD_FWHM_KMS:
            profile = component['amp'] * np.exp(
                -0.5 * ((grid - component['v']) / component['sigma']) ** 2
            )
            # 隐藏峰值 0.5% 以下的尾部，避免虚线长时间贴在基线上。
            visible = np.abs(profile) > 0.005 * abs(component['amp'])
            axes.plot(grid, np.where(visible, profile * multiplier + offset, np.nan),
                      color=color, ls='--', lw=BROAD_COMPONENT_LINEWIDTH)
    axes.plot(grid, model_at(result, grid) * multiplier + offset,
              color=color, ls=':', lw=1.05)


def legend_required_top(figure, axes, legend):
    """计算当前图例位置所需的 y 上限，不改变坐标范围。"""
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    gap = renderer.points_to_pixels(LEGEND_CLEARANCE_POINTS)
    box = legend.get_window_extent(renderer)
    left, right = axes.transData.inverted().transform(
        [(box.x0 - gap, box.y0), (box.x1 + gap, box.y0)]
    )[:, 0]
    left, right = max(left, axes.get_xlim()[0]), min(right, axes.get_xlim()[1])
    peaks = []
    for line in axes.lines:
        # 使用路径顶点，包含阶梯谱线的水平段和竖直段。
        vertices = line.get_path().vertices
        x, y = vertices[:, 0], vertices[:, 1]
        segments = ((np.minimum(x[:-1], x[1:]) <= right)
                    & (np.maximum(x[:-1], x[1:]) >= left))
        values = np.r_[y[:-1][segments], y[1:][segments]]
        values = values[np.isfinite(values)]
        if values.size:
            peaks.append(float(values.max()))
    fraction = (box.y0 - gap - axes.bbox.y0) / axes.bbox.height
    if not peaks or fraction <= 0:
        return axes.get_ylim()[1]
    bottom, top = axes.get_ylim()
    return bottom + (max(peaks) - bottom) / fraction


def reserve_legend_headroom(figure, axes, legend, layout):
    """选择所需留白较少的顶部角落，并按需扩展 y 上限。"""
    locations = ('upper left', 'upper right') if LEGEND_AUTO_SIDE else (LEGEND_LOCATION,)
    candidates = []
    for location in locations:
        legend.set_loc(location)
        candidates.append((legend_required_top(figure, axes, legend), location))
    required_top, location = min(candidates, key=lambda candidate: candidate[0])
    legend.set_loc(location)
    layout['legend_location'] = location
    bottom, top = axes.get_ylim()
    if not AUTO_LEGEND_HEADROOM:
        return
    if required_top > top:
        axes.set_ylim(bottom, required_top)
        layout['ylim'] = [bottom, required_top]


def decorate_profile(figure, axes, source, layout):
    lower, upper = layout['xlim']
    for family_name, offset in layout['baselines_K'].items():
        axes.axhline(offset, color='.65', lw=0.6, zorder=0)
        family_label = family_name.replace('HCO+', 'HCO$^+$')
        if family_name == ARCHIVE_FAMILY and 'transition_labels' in layout:
            family_label = 'archive ' + layout['transition_labels']['12CO10']
        # 基线标签靠右放置，为中部的谱线留出空间。
        axes.text(upper - 0.01 * (upper - lower), offset,
                  family_label + f': {offset:g} K',
                  ha='right', va='center', fontsize=PROFILE_FONT_SIZES['baseline'],
                  bbox=dict(facecolor='white', edgecolor='none', pad=1, alpha=0.9))
    axes.set(xlim=layout['xlim'], ylim=layout['ylim'],
             xlabel=r'$v_{\rm LSR}$ (km s$^{-1}$)', ylabel=r'$T_{\rm mb}$ (K)')
    axes.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=7))
    axes.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=8))
    axes.xaxis.label.set_size(PROFILE_FONT_SIZES['axis_label'])
    axes.yaxis.label.set_size(PROFILE_FONT_SIZES['axis_label'])
    axes.tick_params(direction='in', top=True, right=True,
                     labelsize=PROFILE_FONT_SIZES['tick'])
    axes.set_title(source, fontsize=PROFILE_FONT_SIZES['title'], pad=10)
    handles, labels = axes.get_legend_handles_labels()
    entries = sorted(zip(handles, labels),
                     key=lambda entry: LEGEND_ORDER.index(entry[0].get_gid()))
    handles, labels = zip(*entries)
    legend = axes.legend(handles, labels, loc=LEGEND_LOCATION,
                         ncol=1, fontsize=PROFILE_FONT_SIZES['legend'], frameon=True, facecolor='white',
                         edgecolor='0.35', framealpha=0.95, columnspacing=1, handlelength=2.4)
    for handle in legend.legend_handles:
        handle.set_linewidth(1.5)
    # figure.text(0.5, 0.012,
    #             f'Solid: observed   Dotted: total fit   Dashed: FWHM > {BROAD_FWHM_KMS} km/s',
    #             ha='center', fontsize=PROFILE_FONT_SIZES['caption'])
    figure.subplots_adjust(left=0.14, right=0.97, bottom=0.12, top=0.94)
    reserve_legend_headroom(figure, axes, legend, layout)


def write_profiles(results, sources, layouts, output_directory=None):
    output_directory = Path(output_directory) if output_directory is not None else FIGURE_DIRECTORY
    plt.rcParams.update({'font.family': 'serif', 'font.size': PROFILE_FONT_SIZES['base'],
                         'axes.linewidth': 0.8, 'pdf.fonttype': 42})
    output_directory.mkdir(parents=True, exist_ok=True)
    for source in sources:
        data = sorted([result for result in results if result['source'] == source
                       and result['line'] in layouts[source]['multipliers']],
                      key=lambda result: ORDER.index(result['line']))
        figure, axes = plt.subplots(figsize=FIGURE_SIZE)
        for result in data:
            plot_transition(axes, result, layouts[source])
        decorate_profile(figure, axes, source, layouts[source])
        # figure.savefig(output_directory / f'{source}.pdf')
        figure.savefig(output_directory / f'{source}.png', dpi=PNG_DPI)
        plt.close(figure)


def fmt(value, error):
    """按误差量级确定小数位，让数值和误差的小数位对齐。"""
    if error <= 0:
        return f'${value:.3g}$'
    decimals = min(5, max(0, 1 - math.floor(math.log10(error))))
    return f'${value:.{decimals}f} \\pm {error:.{decimals}f}$'


def write_table(rows):
    heading = r'Target & Transition & $v_0$ & $T_p$ & FWHM & $\int T\,dv$ \\'
    lines = [
        r'\begin{longtable}{cccccc}',
        r'\caption{Fitted broad molecular-line components '
        r'($\mathrm{FWHM}>8$ km s$^{-1}$).}\label{tab:broad}\\',
        r'\toprule', heading,
        r' & & (km s$^{-1}$) & (K) & (km s$^{-1}$) & (K km s$^{-1}$) \\',
        r'\midrule\endfirsthead',
        r'\multicolumn{6}{c}{Table \thetable\ continued}\\',
        r'\toprule ' + heading,
        r'\midrule\endhead',
        r'\midrule\multicolumn{6}{c}{Continued on next page}\\\endfoot',
        r'\bottomrule\endlastfoot',
    ]
    previous = None
    for row_index, row in enumerate(rows):
        if previous and previous != row['source']:
            lines.append(r'\addlinespace')
        source = row['source'] + (r'$^{\dagger}$' if row['legacy'] else '')
        cells = [
            source if previous != row['source'] else '',
            LABEL[row['line']].replace('–', '--'),
            fmt(row['v'], row['ev']), fmt(row['amp'], row['eamp']),
            fmt(row['fwhm'], row['efwhm']),
            fmt(row['integrated'], row['eintegrated']),
        ]
        same_target_next = (row_index + 1 < len(rows) and
                            rows[row_index + 1]['source'] == row['source'])
        lines.append(' & '.join(cells) + (r' \\*' if same_target_next else r' \\'))
        previous = row['source']
    lines.extend([
        r'\end{longtable}',
        r'\noindent\footnotesize Notes. All components passing the strict width threshold '
        'are retained, including absorption, which has negative peak temperature and '
        'integrated intensity. Errors are the reported fit uncertainties. '
        r'$^{\dagger}$Legacy W44 products lack temperature-calibration metadata; '
        'their values are preserved. Some other fits divide an input labelled '
        r'$T_{\rm mb}$ by a beam efficiency; absolute intensities require calibration '
        'verification. Selection by width alone does not establish a shocked-gas origin.',
    ])
    fragment = r'\setlength{\LTcapwidth}{\textwidth}' + '\n' + '\n'.join(lines) + '\n'
    (OUT / 'broad_line_table.tex').write_text(fragment)


def update_profile_document(sources):
    """沿用原脚本的独立文档模板，替换其中的谱图章节。"""
    figure_tex = [
        r'\clearpage', r'\section*{Molecular-line profiles}',
        r'\noindent Velocity windows and vertical limits are adjusted for each target. '
        'Transitions of the same molecule share a baseline; the spacing between molecules '
        'is chosen from their observed and fitted peak-to-trough ranges. CO uses '
        r'$s=1$; HCO$^+$, HCN and SiO use $s=10$. The baseline $b$ is labelled in each '
        'panel. Solid curves show observations, dotted curves total fits, and dashed '
        r'curves individual components with FWHM $>8$ km s$^{-1}$.',
        r'\medskip',
    ]
    for row_index in range(0, len(sources), 2):
        figure_tex.append(r'\noindent')
        for index, source in enumerate(sources[row_index:row_index + 2]):
            if index:
                figure_tex.append(r'\hfill')
            figure_tex.extend([
                r'\begin{minipage}[t]{0.49\textwidth}\centering',
                r'\includegraphics[width=\linewidth]{figures/' + source + '.pdf}',
                r'\end{minipage}',
            ])
        figure_tex.append(r'\par')
        if row_index + 2 < len(sources):
            figure_tex.append(r'\clearpage')
    document = OUT / 'broad_line_table_standalone.tex'
    existing = document.read_text()
    marker = r'\clearpage' + '\n' + r'\section*{Molecular-line profiles}'
    assert marker in existing
    prefix = existing.split(marker, 1)[0]
    prefix = prefix.replace('amsmath,graphicx,tikz', 'amsmath,graphicx')
    prefix = re.sub(r'\\definecolor[^\n]*\n', '', prefix)
    document.write_text(prefix + '\n'.join(figure_tex) + '\n' + r'\end{document}' + '\n')
    return document


def write_summary(results, rows, sources, layouts):
    summary = {
        'fit_products': len(results), 'targets': len(sources),
        'broad_targets': len({row['source'] for row in rows}),
        'broad_components': len(rows),
        'emission': sum(row['type'] > 0 for row in rows),
        'absorption': sum(row['type'] < 0 for row in rows),
        'multipliers': FACTOR, 'target_layouts': layouts,
    }
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2))


def write_readme():
    (OUT / 'README.md').write_text("""# 宽线结果

FWHM > 8 km/s，保留全部宽分量，六列表格居中；负峰温、负积分强度表示吸收。原始分量编号与类型保存在 CSV。

每个目标一张共享 xy 轴的叠加图。同一分子的不同跃迁共用基线，不同分子的间距按当前目标所有观测谱及总拟合的峰谷范围计算，保留间隔并按规则取整。CO 基线为0，下方依次是 HCO+、HCN、SiO。CO 倍率固定1，其他分子固定10。强度恢复公式 T=(y-b)/s。

每个目标独立确定坐标范围：速度覆盖所有拟合分量的中心±3σ，再左右留白并取整；3C391SW 按指定采用 -20 至180 km/s。纵轴包住所有平移后的谱线，上下各留3.5%的余量。图例置于绘图区之外，避免遮挡曲线。具体坐标与偏移记录在 summary.json。

当前 TeX 通过 includegraphics 插入 figures 中的 PDF 图片，每行两张；figures 同时保存 PNG。编译与分享时需保留 figures 相对目录。原表格内容保持不变。旧版图片已由新版覆盖。

selection_manifest.json 记录版本和温标疑点；未重新拟合或重标定，W44 旧结果仍缺少温标说明。

重现：MPLCONFIGDIR=.mplconfig /Users/silence/miniconda3/bin/python code/plotting/build_broad_line_products.py
""")


def validate_layouts(layouts):
    for layout in layouts.values():
        families = list(layout['baselines_K'])
        offsets = layout['baselines_K']
        envelopes = layout['envelopes_K']
        for upper, lower in zip(families, families[1:]):
            assert offsets[lower] + envelopes[lower][1] < offsets[upper] + envelopes[upper][0]


def main():
    OUT.mkdir(exist_ok=True)
    selected, versions = select_fit_files()
    results, rows, warnings = collect_results(selected)
    write_inventory(rows, versions, warnings)

    sources = sorted({result['source'] for result in results})
    layouts = {source: layout_for([result for result in results if result['source'] == source])
               for source in sources}
    write_profiles(results, sources, layouts)

    if '--figures-only' in sys.argv:
        summary = json.loads((OUT / 'summary.json').read_text())
        summary.update(multipliers=FACTOR, target_layouts=layouts)
        (OUT / 'summary.json').write_text(json.dumps(summary, indent=2))
        assert all(float(value) % BASELINE_STEP_K == 0
                   for layout in layouts.values() for value in layout['baselines_K'].values())
        print(f'Updated {len(sources)} figures with adaptive multipliers and separate archive CO baselines.')
        return

    write_table(rows)
    document = update_profile_document(sources)
    write_summary(results, rows, sources, layouts)
    write_readme()
    validate_layouts(layouts)
    print(json.dumps({'targets': len(sources), '3C391SW': layouts['3C391SW'],
                      'tex_bytes': document.stat().st_size}, indent=2))


if __name__ == '__main__':
    main()
