# 3C391SW 单源谱图布局

使用 plot_single_line_profile.py 和同目录 3C391SW_apex_profile.json，按 T_display = Tmb × multiplier + baseline_K 绘制。CO 基线为 0 K，APEX CO(2–1) 为 −10 K；HCO+、HCN、SiO 在之前布局上整体下移 5 K，分别为 −20、−35、−40 K。倍率及横轴 −10–180 km/s 不变，纵轴自动计算。这些偏置仅改变显示，不改变效率校准或窗口积分。

运行时通过 --reference 指定配置，并以 --output-dir 指定 output/apex_adaptation/single_source。输出覆盖 3C391SW.png、3C391SW.pdf 和 3C391SW_layout_used.json；布局记录保存实际使用参数和拟合文件。此次已重新生成以上文件，不修改窗口图和积分 JSON。
