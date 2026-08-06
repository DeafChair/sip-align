# sip-align · SIP 定标与跨夜晚对齐

通用的天文图像定标与对齐工具：

1. **SIP 定标**：给没有畸变模型的 FITS 加上带 SIP3 畸变修正的 WCS。
   星表来自 VizieR 在线服务（UCAC4 / Gaia DR3），免注册、不用本地星表；
2. **跨夜晚对齐**：把另一晚的目标帧对齐到参考帧（模板）的网格上，
   用于测光/叠图/查找瞬变源等需要严格对齐的场景。

核心匹配逻辑复用 MHP 的实战代码（星点检测 → 星表投影 → astroalign
精修 → 互匹配 → SIP3 拟合 → 残差验证），不依赖特定望远镜型号。

## 依赖

- Python 3.10+
- `numpy scipy astropy photutils matplotlib`
- MHP 项目目录（脚本自动加入 `sys.path`）
- 在线模式需要网络

## 快速开始

### 1) 单帧 SIP 定标

```powershell
python sip_calibrate.py image.fts --out solve
```

- 如果 FITS 头里已有 WCS，会直接用它当种子做 SIP3 精修；
- 如果没有 WCS，会按 30° 步长试探方位角（12 个方向 × 正反镜像），
  astroalign 能找回 ±25° 内的偏差，所以只要粗指向没错就能解；
- 如果头里没有像素尺度（XPIXSZ/FOCALLEN），用 `--pixscale` 指定：

```powershell
python sip_calibrate.py image.fts --pixscale 1.72 --pa 180 --out solve
```

输出 `solve/image_sip3.fits`（原图 + SIP3 WCS）和 `report.json`。

### 2) 参考帧与目标帧都定标

分别对模板夜和目标夜的帧跑第 1 步，得到两个 `*_sip3.fits`。

### 3) 跨夜晚对齐

```powershell
python align_nights.py template_sip3.fits target_sip3.fits --out align
```

输出 `align/target_aligned_to_template.fits` 和 `report.json`。
看 `quality_grade`（GOOD）和 `rms_px`（小于 1 px 可直接使用）。

## 参数

### sip_calibrate.py

| 参数 | 说明 |
|---|---|
| `--catalog` | `ucac4`（默认）/ `gaia3`（VizieR 在线）/ `gaia`（本地索引） |
| `--pixscale` | 像素尺度（arcsec/px），头里算不出时必填 |
| `--pa` | 方位角试探列表，默认 0~330 步长 30 |
| `--seed-wcs` | 用另一张已定标帧的 WCS 当种子 |
| `--vizier-mirror` | 手动指定 VizieR 镜像（默认自动切换 CDS/斯特拉斯堡/哈佛） |
| `--gaia-root` | 本地 Gaia 索引目录 |
| `--translation-search` | 星表平移搜索半径（px），指向偏差大时调大 |

### align_nights.py

| 参数 | 说明 |
|---|---|
| `template` | 参考帧（对齐后的网格） |
| `target` | 另一晚的目标帧 |
| `--out` | 输出目录 |

## 测试记录

- 真实巡天帧两遍曝光：定标 SIP3 全场 RMS 0.23 px（在线 UCAC4），
  对齐 2029 颗星 / RMS 0.081 px / GOOD；
- 在线 UCAC4 与本地 Gaia 解全场互差 0.07~0.17 角秒，可直接互换；
- SIP3 效果：视场角点平均误差从 5.15 px 降到 0.25 px（约 20 倍）。

## 常见问题

**为什么用 VizieR？**
公开接口、免注册，Astrometrica 的在线星表用的就是它；astrometry.net
需要 API key 且服务不稳定。

**SIP3 修的是什么？**
修的是“像素 ↔ 天球坐标”的映射（WCS），不修改图片本身；坐标模型修对了，
跨夜晚重投影对齐时边缘才严格对得上。

**头里已经有 WCS 还需要跑吗？**
建议跑。很多初版 WCS 是纯线性模型，没有 SIP 畸变项；`sip_calibrate.py`
会用已有 WCS 当种子，补齐 SIP3 并验证。

**跨夜晚对齐要注意什么？**
两帧都先定标；不同夜晚指向/旋转不同没有关系，WCS 重投影会统一到模板网格，
星表微调再吃掉残余平移。

## 致谢

复用 MHP 项目实战代码；VizieR 为 CDS 公开服务。
