# HCN 窗口谱适配

plot_window_profiles.py 在原有五条谱线之后接受 HCN10、HCN32，沿用现有谱读取、标定、颜色及标签。两条新线的显示偏移分别为 -15、-20 K；参考 JSON 中的倍率优先，缺少 HCN 倍率时使用 1。原有谱线倍率和偏移不变。仅处理实际收集到的谱线，不虚构缺失谱线。

积分仍使用未缩放谱的通道覆盖权重，边缘通道按窗口重叠比例积分；误差仍由 RMS 和权重平方和传播。缺少 RMS 时保留 null 并给出原有警告。HCN 的积分、误差、来源及显示参数写入既有 lines 字典，ratios 仍以 CO43 为分母。输出仍为输出目录中的 {source}_window_profiles.png 和 {source}_window.json，同名文件沿用原覆盖行为。CLI、默认路径和 JSON 层级未改变，W44 目标沿用大写规范。验证产物位于 tmp/hcn_profiles/。
