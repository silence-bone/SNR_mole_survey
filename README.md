# SNR molecular line survey

本仓库当前保存 `code/plotting/` 中的绘图脚本、配置、公共字体设置和已有使用说明，作为后续修改的版本起点。观测 FITS、生成图片、科学结果、临时文件及 Python 缓存不纳入本次版本管理。目录结构与本地项目保持一致，已有科学脚本未因上传而改写。

## FITS 谱线绘图

`plot_fits_spectra.py` 可在同一坐标系叠加多条一维 FITS 谱线，默认速度窗口为 −100 至 200 km/s，支持自定义 RGB 颜色与图例，字体复用 window profiles 的公共设置。参数集中在脚本开头，也可从命令行传入文件。请以当前脚本顶部的参数值为准。

安装 NumPy、Astropy 和 Matplotlib 后，从仓库根目录运行，例如：

```bash
python code/plotting/plot_fits_spectra.py /path/to/line1.fits /path/to/line2.fits \
  --window -100 200 \
  --colors 40,90,210 220,65,45 \
  --labels 'Line 1' 'Line 2'
```

FITS 需要由使用者另外提供。图片默认保存到 `output/fits_spectra/`，使用时间戳命名，然后交给本机系统图片查看器打开；旧图片无需关闭即可再次运行。只保存时加 `--no-show`。详细读取规则、参数与输出行为见 [FITS 谱线绘图说明](code/plotting/plot_fits_spectra.review.md)。

## 其他绘图与配置

`code/plotting/` 还包括单源谱线、window profiles、宽线比值与表格、W44 射电和 H₂ 叠加图等脚本。相关说明保留在同目录的 `.review.md` 中，手动配置位于 `configs/`。公共字体设置位于 `profile_plot_style.py`，另保留已有 `matplotlibrc`。

这次上传是绘图目录的源码快照，尚未包含完整分析流程。部分脚本依赖 `code/other/` 的波束效率和重标定工具、`code/fitting/` 的拟合模块，以及本地数据、历史拟合结果或绝对路径。这些脚本需要完整本地项目及相应输入才能运行，单独克隆当前仓库不足以重现全部结果。除 FITS 绘谱所需依赖外，其他绘图还可能使用 SciPy、Pillow 等，具体见各脚本的导入和已有说明。

## 修改历史与保存范围

GitHub 的提交历史记录首次上传后的每次版本变化，可在文件页面的 History 查看，也可在本地运行：

```bash
git log --oneline
git log --oneline -- code/plotting/plot_fits_spectra.py
git show <commit>
```

目前根目录 `.gitignore` 使用目录白名单，仅开放本说明、忽略规则和 `code/plotting/`，防止把本地的大量数据与结果误加入仓库。以后纳入其他模块时，需要相应扩展白名单。

本次上传仅整理版本管理和仓库说明，未修改绘图算法、参数或历史科学结果。上传前检查 Python 语法与 JSON 格式，上传后按 Git blob 校验源码内容一致性；这不代表所有依赖外部数据的分析流程已经重新运行。
