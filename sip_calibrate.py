#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SIP 定标工具（通用）：给 FITS 加上带 SIP3 畸变模型的 WCS。

原理：
1. 从 FITS 头读指向/像素尺度，或直接使用已有的 WCS / --seed-wcs 当种子；
2. 从 VizieR 在线星表（UCAC4 / Gaia DR3，免注册）下载视场星；
3. 图像找星 → astroalign 精修（方位角有小偏差也能找回）→ 星表匹配；
4. 拟合线性 WCS + SIP3 畸变模型，验证残差后写回 FITS。

用法：
    python sip_calibrate.py image.fts --out out
    python sip_calibrate.py image.fts --pixscale 1.72 --pa 180 --out out
    python sip_calibrate.py image.fts --seed-wcs solved.fits --out out

输出：<图像名>_sip3.fits（原图 + SIP3 WCS）和 report.json。
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parents[1]
MHP_ROOT = ROOT / "MHP"
if str(MHP_ROOT) not in sys.path:
    sys.path.insert(0, str(MHP_ROOT))

from m31_hmtproject.astrometry import calibrate_wcs_with_gaia  # noqa: E402
from m31_hmtproject.grappa3e_catalog import discover_default_catalog_root  # noqa: E402


def parse_hms_dms(text: str) -> float:
    text = str(text).strip()
    parts = text.replace(":", " ").split()
    if len(parts) == 3:
        h, m, s = (float(v) for v in parts)
        return (h + m / 60.0 + s / 3600.0) * 15.0
    return float(text)


def parse_dms(text: str) -> float:
    text = str(text).strip()
    sign = -1.0 if text.startswith("-") else 1.0
    text = text.lstrip("+-")
    parts = text.replace(":", " ").split()
    if len(parts) == 3:
        d, m, s = (float(v) for v in parts)
        return sign * (d + m / 60.0 + s / 3600.0)
    return sign * float(text)


def header_center(header) -> tuple[float, float]:
    """RA/Dec of the pointing in degrees from common header keywords."""
    ra = None
    dec = None
    for ra_key in ("OBJCTRA", "RA"):
        if ra_key in header and header[ra_key]:
            ra = parse_hms_dms(header[ra_key])
            break
    for dec_key in ("OBJCTDEC", "DEC"):
        if dec_key in header and header[dec_key]:
            dec = parse_dms(header[dec_key])
            break
    if ra is None or dec is None:
        raise ValueError("FITS 头里找不到 RA/DEC（OBJCTRA/OBJCTDEC 或 RA/DEC）")
    return ra, dec


def auto_pixel_scale(header) -> float | None:
    """arcsec/px from XPIXSZ(um) and FOCALLEN(mm)."""
    try:
        xpix = float(header.get("XPIXSZ", 0.0) or 0.0)
        focal = float(header.get("FOCALLEN", 0.0) or 0.0)
        if xpix > 0 and focal > 0:
            return 206265.0 * (xpix * 1e-6) / (focal * 1e-3)
    except (TypeError, ValueError):
        pass
    return None


def header_has_wcs(header) -> bool:
    has_basic = all(k in header for k in ("CTYPE1", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2"))
    has_matrix = ("CD1_1" in header) or (("PC1_1" in header) and ("CDELT1" in header))
    return has_basic and has_matrix


def build_seed_header(header, ra_deg: float, dec_deg: float, scale_arcsec: float, pa_deg: float, flip: bool):
    seed = header.copy()
    w, h = int(seed.get("NAXIS1", 4096)), int(seed.get("NAXIS2", 4096))
    s = scale_arcsec / 3600.0
    a = np.radians(pa_deg)
    cos_a, sin_a = np.cos(a), np.sin(a)
    c11 = -s * cos_a
    c12 = s * sin_a
    c21 = s * sin_a
    c22 = s * cos_a
    if flip:
        c11, c12 = -c11, -c12
    seed["CTYPE1"] = "RA---TAN"
    seed["CTYPE2"] = "DEC--TAN"
    seed["CUNIT1"] = "deg"
    seed["CUNIT2"] = "deg"
    seed["CRPIX1"] = w / 2.0 + 0.5
    seed["CRPIX2"] = h / 2.0 + 0.5
    seed["CRVAL1"] = ra_deg
    seed["CRVAL2"] = dec_deg
    seed["CD1_1"] = c11
    seed["CD1_2"] = c12
    seed["CD2_1"] = c21
    seed["CD2_2"] = c22
    seed["RADESYS"] = "ICRS"
    seed["EQUINOX"] = 2000.0
    return seed


def safe_header(d: dict, base=None) -> fits.Header:
    """Build a FITS Header from a dict, skipping HISTORY/COMMENT entries."""
    out = fits.Header(base) if base is not None else fits.Header()
    for key, value in d.items():
        if not key or key in {"HISTORY", "COMMENT", ""}:
            continue
        try:
            out[key] = value
        except Exception:
            pass
    return out


def trial_result_pretty(result) -> dict:
    return {
        "status": result.status,
        "selected_model": result.selected_model,
        "matched_sources": result.matched_sources,
        "match_rms_px": result.match_rms_px,
        "sip_rms_px": result.sip_rms_px,
        "sip_validation_accepted": result.sip_validation_accepted,
        "sip_edge_rms_px": result.sip_edge_rms_px,
        "sip_corner_rms_px": result.sip_corner_rms_px,
        "sip_rejection_reasons": result.sip_rejection_reasons or [],
        "notes": result.notes,
    }


def score_trial(result) -> tuple[float, int]:
    """Lower is better; ok+sip beats ok+linear beats rejected."""
    if result.status not in {"ok", "qc_ok"}:
        return (100.0, 0)
    model_rank = 0 if str(result.selected_model).startswith("sip") else 1
    rms = result.sip_rms_px if result.sip_rms_px is not None else result.match_rms_px
    if rms is None:
        rms = 999.0
    return (model_rank, -float(rms))


def solve_one(
    data: np.ndarray,
    header,
    catalog_root: Path,
    out_dir: Path,
    *,
    pa_list,
    flips,
    mode: str,
    sip_degree: int,
    refine_prior: bool = True,
    pixscale_override: float | None = None,
):
    ra_deg, dec_deg = header_center(header)
    scale = pixscale_override if pixscale_override is not None else auto_pixel_scale(header)
    if scale is None:
        raise RuntimeError("无法从 FITS 头算像素尺度，请用 --pixscale 指定（arcsec/px）")
    print(f"[SEED] RA={ra_deg:.6f} Dec={dec_deg:+.6f} scale={scale:.4f}\"/px "
          f"PA={pa_list} flip={flips}")
    best = None
    trials = []
    refined_best = None
    for pa in pa_list:
        for flip in flips:
            seed = build_seed_header(header, ra_deg, dec_deg, scale, pa, flip)
            tag = f"pa{int(pa):03d}{'_flip' if flip else ''}"
            if refine_prior:
                from m31_hmtproject.next_astrometry import _refine_camera_prior_with_gaia

                try:
                    refined, matched, rms, pixscale, parity = _refine_camera_prior_with_gaia(
                        data, seed, str(catalog_root)
                    )
                except Exception as exc:
                    trials.append({
                        "pa": pa, "flip": flip, "refine": "failed",
                        "exception": f"{type(exc).__name__}: {exc}",
                    })
                    print(f"[REFINE] PA={pa:6.1f} flip={flip} FAIL {type(exc).__name__}: {exc}")
                    continue
                trials.append({
                    "pa": pa, "flip": flip, "refine": "ok",
                    "matched": int(matched), "rms_px": float(rms),
                    "pixel_scale": float(pixscale), "parity": parity,
                })
                print(f"[REFINE] PA={pa:6.1f} flip={flip} OK matched={matched} rms={rms:.3f}px")
                if refined_best is None or (int(matched), -float(rms)) > (
                    refined_best[1], -refined_best[2]
                ):
                    refined_best = (tag, int(matched), float(rms), refined)
                continue
            trial_dir = out_dir / tag
            trial_dir.mkdir(parents=True, exist_ok=True)
            try:
                result = calibrate_wcs_with_gaia(
                    data, seed, str(catalog_root), trial_dir,
                    mode=mode, sip_degree=sip_degree,
                )
            except Exception as exc:
                trials.append({"pa": pa, "flip": flip, "exception": f"{type(exc).__name__}: {exc}"})
                print(f"[TRIAL] PA={pa:6.1f} flip={flip} EXC {type(exc).__name__}: {exc}")
                continue
            pretty = trial_result_pretty(result)
            pretty["pa"] = pa
            pretty["flip"] = flip
            trials.append(pretty)
            print(f"[TRIAL] PA={pa:6.1f} flip={flip} status={result.status} "
                  f"model={result.selected_model} match={result.match_rms_px} "
                  f"sip={result.sip_rms_px}")
            if result.status in {"ok", "qc_ok"} and (
                best is None or score_trial(result) < score_trial(best[1])
            ):
                best = (trial_dir, result)
    if refine_prior and refined_best is not None:
        tag, matched, seed_rms, refined_header = refined_best
        trial_dir = out_dir / "refine_best" / tag
        trial_dir.mkdir(parents=True, exist_ok=True)
        try:
            result = calibrate_wcs_with_gaia(
                data, refined_header, str(catalog_root), trial_dir,
                mode=mode, sip_degree=sip_degree,
            )
        except Exception as exc:
            trials.append({"refined_seed": tag, "exception": f"{type(exc).__name__}: {exc}"})
            print(f"[SIP3] {tag} EXC {type(exc).__name__}: {exc}")
            return None, trials
        pretty = trial_result_pretty(result)
        pretty["refined_seed"] = tag
        pretty["refined_matched"] = matched
        pretty["refined_rms_px"] = seed_rms
        trials.append(pretty)
        print(f"[SIP3] {tag} status={result.status} model={result.selected_model} "
              f"match={result.match_rms_px} sip={result.sip_rms_px}")
        if result.status in {"ok", "qc_ok"}:
            best = (trial_dir, result)
    return best, trials


def main() -> int:
    parser = argparse.ArgumentParser(description="Generic SIP3 WCS calibration (VizieR online / local Gaia)")
    parser.add_argument("fits", help="FITS image to calibrate")
    parser.add_argument("--out", default="sip_calibrate_output", help="output directory")
    parser.add_argument("--catalog", default="ucac4", choices=["ucac4", "gaia3", "gaia"],
                        help="ucac4/gaia3 = VizieR online (default ucac4); gaia = local index")
    parser.add_argument("--gaia-root", default=None, help="local compact Gaia EDR3 index root")
    parser.add_argument("--vizier-mirror", default=None,
                        help="override VizieR mirror base URL, e.g. https://vizier.cfa.harvard.edu")
    parser.add_argument("--mode", default="sip3", choices=["qc", "linear", "sip3"])
    parser.add_argument("--sip-degree", type=int, default=3)
    parser.add_argument("--pixscale", type=float, default=None,
                        help="pixel scale in arcsec/px (auto from XPIXSZ/FOCALLEN if possible)")
    parser.add_argument("--pa", default="0,30,60,90,120,150,180,210,240,270,300,330",
                        help="position-angle trials in degrees (only used when no WCS/seed)")
    parser.add_argument("--seed-wcs", default=None,
                        help="use WCS from another solved FITS as the starting point")
    parser.add_argument("--no-flip", action="store_true", help="disable mirrored seed trials")
    parser.add_argument("--translation-search", type=float, default=40.0,
                        help="star-table translation search radius in px")
    return _run(parser.parse_args())


def _run(args) -> int:
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.catalog != "gaia":
        from vizier_catalog import VizieRCatalog
        from m31_hmtproject import astrometry as mhp_astrometry
        from m31_hmtproject import next_astrometry as mhp_next_astrometry

        mirrors = [args.vizier_mirror] if args.vizier_mirror else None
        factory = lambda root: VizieRCatalog(root, source=args.catalog, mirrors=mirrors)
        mhp_astrometry.Grappa3ECatalog = factory
        mhp_next_astrometry.Grappa3ECatalog = factory
        catalog_root = out_dir / "vizier-online"
        print(f"[CATALOG] 在线 {args.catalog}（VizieR ASU，不需要注册/本地星表）")
    else:
        catalog_root = Path(args.gaia_root) if args.gaia_root else discover_default_catalog_root()
        if catalog_root is None or not catalog_root.is_dir():
            raise RuntimeError("找不到本地 Gaia 目录；用 --gaia-root 指定，或 --catalog ucac4 使用在线星表")
        print(f"[CATALOG] 本地 Gaia 索引: {catalog_root}")

    fits_path = Path(args.fits)
    with fits.open(fits_path, memmap=False) as hdul:
        data = np.asarray(hdul[0].data, dtype=np.float32)
        header = hdul[0].header

    seed_header = None
    seed_source = None
    if args.seed_wcs:
        with fits.open(args.seed_wcs, memmap=False) as hdul:
            seed_header = hdul[0].header
        seed_source = f"seed-wcs:{args.seed_wcs}"
    elif header_has_wcs(header):
        seed_header = header.copy()
        seed_source = "header-wcs"

    if seed_header is not None:
        work = seed_header.copy()
        for key in ("NAXIS1", "NAXIS2"):
            if key in work:
                work.pop(key, None)
        for key in ("NAXIS1", "NAXIS2"):
            work[key] = header.get(key, 0)
        trial_dir = out_dir / "solve" / "seed"
        trial_dir.mkdir(parents=True, exist_ok=True)
        result = calibrate_wcs_with_gaia(
            data, work, str(catalog_root), trial_dir,
            mode=args.mode, sip_degree=args.sip_degree,
            translation_search_px=args.translation_search,
        )
        trials = [trial_result_pretty(result)]
        best = (trial_dir, result)
        print(f"[SEED] {seed_source}: {result.status} model={result.selected_model} "
              f"match={result.match_rms_px} sip={result.sip_rms_px}")
    else:
        pa_list = [float(v) for v in str(args.pa).split(",") if v.strip()]
        flips = [False] if args.no_flip else [False, True]
        best, trials = solve_one(
            data, header, catalog_root, out_dir / "solve",
            pa_list=pa_list, flips=flips,
            mode=args.mode, sip_degree=args.sip_degree,
            refine_prior=True,
            pixscale_override=args.pixscale,
        )

    if best is None:
        raise RuntimeError("所有定标尝试都失败；看输出目录里的 trial 日志")
    trial_dir, result = best
    final_header = safe_header(result.header, base=fits.Header(header))
    first_out = out_dir / f"{fits_path.stem}_sip3.fits"
    fits.writeto(first_out, data, final_header, overwrite=True)
    print(f"[OUT] {first_out} | model={result.selected_model} sip_rms={result.sip_rms_px}px")

    report = {
        "input": str(fits_path),
        "seed_source": seed_source,
        "catalog": args.catalog,
        "output": str(first_out),
        "best": trial_result_pretty(result),
        "trials": trials,
    }
    (out_dir / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[REPORT] {out_dir / 'report.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
