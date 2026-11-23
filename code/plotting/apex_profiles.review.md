# APEX CO(2–1) 绘图适配

`build_broad_line_products.py` 将 `_apex_fit.txt` 与普通拟合独立分组选版，内部跃迁键为 `12CO21_apex`，不会覆盖原有 `12CO21`。原始谱补充读取也识别 `_apex.fits` 和 `_apex_re.fits`。拟合元数据中的输入路径、效率、合并通道和基线处理沿用原流程；无拟合时由 FITS 温标和波束效率读取器决定校准。显示名称为 `12CO(2-1) (apex)`，配独立颜色。

单源 profile 的参考 JSON 在 `multipliers` 中设置 `12CO21_apex` 即可纳入绘制；该观测使用独立 `apex CO(2–1)` 基线，缺少基线设置时自动补齐，与 archive 分层方式一致。批量自动布局也将该观测单独放置，倍率为 1。窗口 profile 默认纳入这一观测，缺少数据时跳过，与其他跃迁相同。窗口积分 JSON 以独立键保存积分、误差、倍率、偏置和输入元数据，不会混合不同望远镜的 CO(2–1)。所有温度偏置仅用于展示，积分仍使用未乘展示倍率的谱。

输出沿用现有命名：单源目录中的 `<source>.png`、`<source>.pdf` 和 `<source>_layout_used.json`；窗口目录中的 `<source>_window_profiles.png` 和 `<source>_window.json`，同名输出覆盖，W44 区域名保持大写。测试产物放在 `tmp/apex_plot_check/`。已实际运行 3C391SW 的窗口绘图和带 APEX 配置的单源绘图。
