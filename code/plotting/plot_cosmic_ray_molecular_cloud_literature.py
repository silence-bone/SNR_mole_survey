from __future__ import annotations

import csv
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "cosmic_rays_molecular_clouds_representative_sample.csv"
OUTPUT = ROOT / "cosmic_rays_molecular_clouds_timeline.png"
FONT_PATH = "/System/Library/Fonts/Hiragino Sans GB.ttc"

CHANNEL_ORDER = [
    "Microphysics", "Thermal", "Dust_chemistry", "UV_chemistry", "Transport",
    "Acceleration", "mm_submm", "FIR_submm", "NIR_chemistry", "NIR_H2",
    "Xray", "Gamma", "Thermal_chemistry", "Multiwavelength",
]

CHANNEL_LABELS = {
    "Microphysics": "电离微物理", "Thermal": "热平衡", "Dust_chemistry": "尘埃/冰幔",
    "UV_chemistry": "二次 UV/化学", "Transport": "输运/屏蔽", "Acceleration": "局地加速",
    "mm_submm": "毫米/亚毫米离子", "FIR_submm": "远红外/亚毫米离子",
    "NIR_chemistry": "近红外 H3+", "NIR_H2": "近红外 CR-H2",
    "Xray": "X 射线 Fe Kα", "Gamma": "γ 射线", "Thermal_chemistry": "热化学/CRDR",
    "Multiwavelength": "多波段综述",
}

ENV_COLORS = {
    "Nearby_clouds": "#2a9d8f", "Quiescent_cloud": "#2a9d8f", "Cold_cloud": "#2a9d8f",
    "Starless_core": "#2a9d8f", "Dense_core": "#457b9d", "Dense_cloud": "#457b9d",
    "Dense_clumps": "#457b9d", "High_mass_clumps": "#457b9d", "Diffuse_cloud": "#8ab17d",
    "Diffuse_ISM": "#8ab17d", "Galactic_disk": "#8ab17d", "Galactic_center": "#e76f51",
    "Galactic_center_cluster": "#e76f51", "Star_forming_complex": "#f4a261",
    "Protostar": "#f4a261", "Protostellar_shock": "#f4a261", "Starburst": "#9b5de5",
    "Star_forming_galaxies": "#9b5de5", "General_clouds": "#6c757d", "General": "#6c757d",
}

ANNOTATE = {
    "Spitzer_Tomasko", "Prasad_Tarafdar", "Issa_Wolfendale", "HESS_GC_ridge",
    "Indriolo_McCall", "Fermi_Orion", "HESS_PeVatron", "HAWC_Cygnus",
    "Yang_et_al_shielding", "Bialy_et_al_B68", "Neufeld_et_al_B68",
}

ANNOTATE_OFFSETS = {
    "Spitzer_Tomasko": (16, -48),
    "Prasad_Tarafdar": (16, -48),
    "Issa_Wolfendale": (16, -45),
    "HESS_GC_ridge": (16, -48),
    "Indriolo_McCall": (16, -58),
    "Fermi_Orion": (16, -62),
    "HESS_PeVatron": (16, -62),
    "HAWC_Cygnus": (-20, -64),
    "Yang_et_al_shielding": (-190, -20),
    "Bialy_et_al_B68": (-235, -58),
    "Neufeld_et_al_B68": (-225, -24),
}


def font(size: int, index: int = 0) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(FONT_PATH, size=size, index=index)


def polygon_star(cx: float, cy: float, outer: float, inner: float, points: int = 5):
    coords = []
    for i in range(points * 2):
        angle = -math.pi / 2 + i * math.pi / points
        radius = outer if i % 2 == 0 else inner
        coords.append((cx + radius * math.cos(angle), cy + radius * math.sin(angle)))
    return coords


def draw_marker(draw: ImageDraw.ImageDraw, x: float, y: float, marker: str, color: str, r: int = 12) -> None:
    outline = "#ffffff"
    width = 3
    if marker == "o":
        draw.ellipse((x-r, y-r, x+r, y+r), fill=color, outline=outline, width=width)
    elif marker == "D":
        draw.polygon([(x, y-r-2), (x+r+2, y), (x, y+r+2), (x-r-2, y)], fill=color, outline=outline)
    elif marker == "^":
        draw.polygon([(x, y-r-3), (x+r+3, y+r+2), (x-r-3, y+r+2)], fill=color, outline=outline)
    elif marker == "s":
        draw.rectangle((x-r, y-r, x+r, y+r), fill=color, outline=outline, width=width)
    else:
        draw.polygon(polygon_star(x, y, r+6, (r+6)*0.45), fill=color, outline=outline)


def marker_for(work_type: str) -> str:
    return {
        "Observation": "o", "Observation_model": "D", "Model_observation": "D",
        "Theory": "^", "Model": "^", "Simulation": "s", "Review": "*",
    }.get(work_type, "o")


def main() -> None:
    with INPUT.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    available = [c for c in CHANNEL_ORDER if any(r["channel"] == c for r in rows)]
    width, height = 3400, 2200
    left, right, top, bottom = 510, 120, 300, 430
    plot_w = width - left - right
    plot_h = height - top - bottom
    lane_h = plot_h / max(1, len(available) - 1)
    y_for = {channel: top + idx * lane_h for idx, channel in enumerate(available)}

    def x_for(year: int) -> float:
        return left + (year - 1965) / (2028 - 1965) * plot_w

    image = Image.new("RGB", (width, height), "#ffffff")
    draw = ImageDraw.Draw(image)
    draw.rounded_rectangle((35, 35, width-35, height-35), radius=28, fill="#fbfaf7", outline="#e7e3da", width=2)

    draw.text((left, 80), "宇宙线—分子云研究的跨波段演进（代表性样本）", font=font(55), fill="#172026")
    draw.text(
        (left, 160),
        "点色表示环境；点形表示工作类型。该图展示研究窗口的演进，不代表完整 ADS 文献数量。",
        font=font(27), fill="#5f6368",
    )

    label_font = font(28)
    tick_font = font(23)
    note_font = font(20)

    for idx, channel in enumerate(available):
        y = y_for[channel]
        draw.line((left, y, width-right, y), fill="#ebe7de", width=2)
        label = CHANNEL_LABELS[channel]
        bbox = draw.textbbox((0, 0), label, font=label_font)
        draw.text((left - 32 - (bbox[2]-bbox[0]), y - 18), label, font=label_font, fill="#29343a")

    for year in range(1970, 2030, 5):
        x = x_for(year)
        draw.line((x, top-15, x, top+plot_h+25), fill="#d9d9d9", width=2)
        text = str(year)
        bbox = draw.textbbox((0, 0), text, font=tick_font)
        draw.text((x-(bbox[2]-bbox[0])/2, top+plot_h+35), text, font=tick_font, fill="#5f6368")

    per_lane_counts: dict[tuple[str, int], int] = {}
    for row in rows:
        year = int(row["year"])
        x = x_for(year)
        key = (row["channel"], year)
        count = per_lane_counts.get(key, 0)
        per_lane_counts[key] = count + 1
        y = y_for[row["channel"]] + (count % 3 - 1) * 17
        color = ENV_COLORS.get(row["environment"], "#6c757d")
        marker = marker_for(row["work_type"])
        draw_marker(draw, x, y, marker, color, r=13)
        if row["short_ref"] in ANNOTATE:
            label = row["short_ref"].replace("_", " ")
            dx, dy = ANNOTATE_OFFSETS.get(row["short_ref"], (16, -45))
            draw.text((x+dx, y+dy), label, font=note_font, fill="#263238", stroke_width=2, stroke_fill="#fbfaf7")

    draw.text((left + plot_w/2 - 70, top+plot_h+92), "发表年份", font=font(27), fill="#36454f")

    legend_y = height - 285
    legend_font = font(23)
    env_legend = [
        ("#2a9d8f", "本地/静态/无星云"), ("#457b9d", "致密核/团块"),
        ("#8ab17d", "弥散云/银盘"), ("#e76f51", "银河中心"),
        ("#f4a261", "恒星形成区/原恒星"), ("#9b5de5", "星暴/外星系"),
        ("#6c757d", "通用/综述"),
    ]
    x = left
    for color, label in env_legend:
        draw_marker(draw, x, legend_y, "o", color, r=11)
        draw.text((x+22, legend_y-16), label, font=legend_font, fill="#3c4448")
        x += 350 if len(label) > 8 else 270

    type_y = height - 190
    type_legend = [("o", "观测"), ("D", "观测+模型"), ("^", "理论/模型"), ("s", "数值模拟"), ("*", "综述")]
    x = left
    for marker, label in type_legend:
        draw_marker(draw, x, type_y, marker, "#5b6570", r=11)
        draw.text((x+22, type_y-16), label, font=legend_font, fill="#3c4448")
        x += 250

    draw.text((width-right-930, height-190), "数据：手工筛选的代表性论文（1968–2026）", font=font(22), fill="#777d80")
    image.save(OUTPUT, format="PNG", optimize=True)
    print(OUTPUT)


if __name__ == "__main__":
    main()
