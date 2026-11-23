"""叠加绘制一维 FITS 谱线，保存后非阻塞地打开图片。"""

import argparse
from datetime import datetime
from pathlib import Path
import re
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')  # 图片查看器负责显示，Python 不等待绘图窗口。
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.constants import c
from astropy.io import fits

try:
    from .profile_plot_style import FONT_SIZES, plot_context, style_axes
except ImportError:
    from profile_plot_style import FONT_SIZES, plot_context, style_axes


# ------------------------- 可直接修改的参数 -------------------------
ROOT = Path(__file__).resolve().parents[2]
# 此处相对路径以项目根目录为基准；命令行输入路径以运行目录为基准。
# 每条谱线各写一个字典；label=None 使用完整文件名，rgb=None 自动选色。
SPECTRA = [
    dict(file='snr_data/3C391SW_12CO21.fits', rgb=(40, 90, 210), label=None),
    dict(file='snr_data/3C391SW_12CO43.fits', rgb=(220, 65, 45), label=None),
]
VELOCITY_WINDOW = (-100, 200)  # km/s
YLIM = None                   # None 自动按窗口内的强度定范围，或 (下限, 上限)
Y_LABEL = None                # None 使用 FITS BUNIT；不同标度会分别注明于图例
TITLE = ''
FIGSIZE = (10, 6)
DPI = 200
LINE_WIDTH = 1.2
LEGEND_LOCATION = 'best'
OUTPUT_DIR = ROOT / 'output/fits_spectra'
OUTPUT_FILE = None            # None 每次用时间戳生成 PNG，避免覆盖
SHOW_FIGURE = True            # 用系统图片查看器打开，打开后可立即再运行
# -------------------------------------------------------------------


def normalizeName(text):
    """输出中的 W44 分区名统一大写，保留实际输入文件路径。"""
    return re.sub(r'w44[wncs]', lambda match: match.group().upper(), str(text),
                  flags=re.IGNORECASE)


def rgbColor(rgb):
    """RGB 使用 0--255 数值，包括小数；不混用 0--1 标度。"""
    if rgb is None:
        return None
    values = np.asarray(rgb, dtype=float)
    if values.shape != (3,) or not np.all(np.isfinite(values)) or np.any(
            (values < 0) | (values > 255)):
        raise ValueError('RGB 必须包含三个 0--255 范围内的数值。')
    return tuple(values / 255)


def velocityAxis(header, nchan):
    """读取 CLASS 速度定义或标准线性速度/频率轴，返回 km/s。"""
    pixels = np.arange(nchan) + 1
    deltaPixels = pixels - header['CRPIX1']
    # 项目 CLASS 导出的频率偏移轴，速度信息另存在这两个卡片中。
    if 'DELTAV' in header and 'VELO-LSR' in header:
        return (header['VELO-LSR'] + deltaPixels * header['DELTAV']) / 1000

    axisType = str(header.get('CTYPE1', '')).strip().upper()
    coordinate = header['CRVAL1'] + deltaPixels * header['CDELT1']
    if axisType in ('VRAD', 'VELO', 'VELO-LSR', 'VLSR', 'VOPT'):
        return (coordinate * u.Unit(header.get('CUNIT1', 'm/s'))).to_value(u.km/u.s)
    if axisType != 'FREQ':
        raise ValueError(f'不支持的谱轴 {axisType!r}；需要线性速度轴或 FREQ 轴。')

    frequency = (coordinate * u.Unit(header.get('CUNIT1', 'Hz'))).to_value(u.Hz)
    restFrequency = header.get('RESTFRQ', header.get('RESTFREQ'))  # Hz
    if restFrequency is None or restFrequency <= 0:
        raise ValueError('频率轴转速度需要正的 RESTFRQ 或 RESTFREQ（Hz）。')
    if header['CRVAL1'] == 0 and 'VELO-LSR' in header:
        return header['VELO-LSR']/1000 - c.to_value(u.km/u.s)*frequency/restFrequency
    return c.to_value(u.km/u.s) * (1 - frequency/restFrequency)  # radio 定义


def readSpectrum(entry):
    """读取主 HDU 的一维谱；去除单元素轴，不对数据立方体做隐式平均。"""
    path = Path(entry['file']).expanduser()
    with fits.open(path, memmap=False) as hdul:
        header = hdul[0].header.copy()
        intensity = np.asarray(hdul[0].data, dtype=float).squeeze()
    if intensity.ndim != 1 or intensity.size != header['NAXIS1']:
        raise ValueError(f'{path.name}: 需要第 1 轴为谱轴的一维谱，实际形状 {intensity.shape}。')
    velocity = velocityAxis(header, intensity.size)
    order = np.argsort(velocity)
    label = path.name if entry.get('label') is None else entry['label']
    return dict(velocity=velocity[order], intensity=intensity[order],
                label=normalizeName(label), color=rgbColor(entry.get('rgb')),
                unit=str(header.get('BUNIT', '')).strip())


@plot_context()
def plotSpectra(spectra, window, outputFile):
    """按输入顺序叠加原始通道谱，保留 NaN 缺口。"""
    fig, ax = plt.subplots(figsize=FIGSIZE)
    units = {row['unit'] for row in spectra}
    visibleValues = []
    for row in spectra:
        velocity, intensity = row['velocity'], row['intensity']
        visible = (velocity >= window[0]) & (velocity <= window[1])
        finite = visible & np.isfinite(intensity) & np.isfinite(velocity)
        if not np.any(finite):
            plt.close(fig)
            raise ValueError(f"{row['label']}: 速度窗口内没有有效通道。")
        visibleValues.append(intensity[finite])
        label = row['label']
        if len(units) > 1:
            label += f" [{row['unit'] or 'unit unspecified'}]"
        ax.step(velocity, intensity, where='mid', lw=LINE_WIDTH,
                color=row['color'], label=label)

    if YLIM is None:
        values = np.concatenate(visibleValues)
        low, high = values.min(), values.max()
        margin = 0.05 * (high-low) if high > low else max(abs(low)*0.05, 0.1)
        ax.set_ylim(low-margin, high+margin)
    else:
        ax.set_ylim(*YLIM)
    unit = spectra[0]['unit'] if len(units) == 1 else ''
    ylabel = Y_LABEL or (f'Intensity ({unit})' if unit else 'Intensity')
    ax.set(xlim=window, xlabel=r'Velocity (km s$^{-1}$)', ylabel=ylabel)
    style_axes(ax, normalizeName(TITLE) if TITLE else None)
    ax.legend(loc=LEGEND_LOCATION, fontsize=FONT_SIZES['spectrum_label'], frameon=False)
    fig.tight_layout()
    outputFile.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outputFile, dpi=DPI)
    plt.close(fig)


def openFigure(path):
    """启动外部查看器，不等待其退出；旧图片无需关闭。"""
    if sys.platform == 'win32':
        import os
        os.startfile(str(path))
    else:
        viewer = 'open' if sys.platform == 'darwin' else 'xdg-open'
        subprocess.Popen([viewer, str(path)], start_new_session=True)


def parseArguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('files', nargs='*', type=Path, help='FITS 文件，按此顺序绘制')
    parser.add_argument('--window', nargs=2, type=float, default=VELOCITY_WINDOW,
                        metavar=('VMIN', 'VMAX'), help='速度窗口，km/s，默认 -100 200')
    parser.add_argument('--colors', nargs='+', help='每条谱的 RGB，例如 40,90,210 220,65,45')
    parser.add_argument('--labels', nargs='+', help='每条谱的图例，含空格时加引号')
    parser.add_argument('--output', type=Path, default=OUTPUT_FILE, help='输出图片路径（同名覆盖）')
    parser.add_argument('--no-show', action='store_true', help='仅保存，不打开图片')
    args = parser.parse_args()
    if not np.all(np.isfinite(args.window)) or args.window[0] >= args.window[1]:
        parser.error('速度窗口需要有限且递增的上下限。')
    entries = [dict(file=path) for path in args.files] if args.files else [
        dict(entry, file=ROOT / Path(entry['file']).expanduser()) for entry in SPECTRA]
    if not entries:
        parser.error('请填写 SPECTRA 或在命令行提供 FITS 文件。')
    for option, key in [('colors', 'rgb'), ('labels', 'label')]:
        values = getattr(args, option)
        if values is not None:
            if len(values) != len(entries):
                parser.error(f'--{option} 的数量必须与 FITS 文件数量一致。')
            for entry, value in zip(entries, values):
                if key == 'rgb':
                    try:
                        value = [float(part) for part in value.split(',')]
                        rgbColor(value)
                    except ValueError as error:
                        parser.error(f'无效 RGB {value!r}: {error}')
                entry[key] = value
    return args, entries


def main():
    args, entries = parseArguments()
    spectra = [readSpectrum(entry) for entry in entries]
    if args.output is None:
        stem = normalizeName(Path(entries[0]['file']).stem)
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        outputFile = OUTPUT_DIR / f'{stem}_spectra_{stamp}.png'
    else:
        path = Path(args.output).expanduser()
        outputFile = path.with_name(normalizeName(path.name))
    outputFile = outputFile.resolve()
    plotSpectra(spectra, args.window, outputFile)
    print(f'Saved: {outputFile}')
    if SHOW_FIGURE and not args.no_show:
        openFigure(outputFile)


if __name__ == '__main__':
    main()
