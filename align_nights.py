#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""跨夜晚对齐：把某一晚的目标帧对齐到参考帧（模板）的网格上。

用法:
    python align_nights.py template.fits target.fits --out align

两帧都建议先经过 sip_calibrate.py 定标（带 SIP3 WCS）。本脚本调用
MHP 的 align_ref_to_new：
    WCS 重投影（含 SIP 畸变修复）→ 星表微调 → 亚像素精修 → 质量门禁
输出目标帧对齐到模板网格的 FITS + report.json。
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
MHP_ROOT = ROOT / "MHP"
if str(MHP_ROOT) not in sys.path:
    sys.path.insert(0, str(MHP_ROOT))

from m31_hmtproject.align import align_ref_to_new  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="Cross-night alignment via MHP align_ref_to_new")
    parser.add_argument("template", help="reference FITS (target grid, e.g. a template night)")
    parser.add_argument("target", help="FITS from another night to align onto the template grid")
    parser.add_argument("--out", default="align_output", help="output directory")
    return _run(parser.parse_args())


def _run(args) -> int:
    template_path = Path(args.template)
    target_path = Path(args.target)
    with fits.open(template_path, memmap=False) as hdul:
        template_data = np.asarray(hdul[0].data, dtype=np.float32)
        template_header = hdul[0].header
    with fits.open(target_path, memmap=False) as hdul:
        target_data = np.asarray(hdul[0].data, dtype=np.float32)
        target_header = hdul[0].header

    result = align_ref_to_new(
        ref=target_data,
        new=template_data,
        ref_header=target_header,
        new_header=template_header,
        prefer_wcs=True,
    )

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    aligned_path = out_dir / f"{target_path.stem}_aligned_to_{template_path.stem}.fits"
    fits.writeto(aligned_path, result.ref_aligned, template_header, overwrite=True)

    report = {
        "template": str(template_path),
        "target": str(target_path),
        "aligned": str(aligned_path),
        "method": result.method,
        "matched_stars": result.matched_stars,
        "rms_px": result.rms_px,
        "valid_fraction": result.valid_fraction,
        "shift_yx": list(result.shift_yx),
        "alignment_ok": result.alignment_ok,
        "quality_grade": result.quality_grade,
        "notes": result.notes,
    }
    (out_dir / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"[ALIGN] method={result.method} matched={result.matched_stars} "
          f"rms={result.rms_px}px valid={result.valid_fraction:.3f} "
          f"grade={result.quality_grade} ok={result.alignment_ok}")
    print(f"[OUT] {aligned_path}")
    print(f"[REPORT] {out_dir / 'report.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
