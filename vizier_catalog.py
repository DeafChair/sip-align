#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""VizieR 在线星表适配器 —— 给 MHP 定标器当 Grappa3ECatalog 的替身。

原理：VizieR 是 CDS 的公开星表服务，Astrometrica 的在线星表用的就是它的
ASU (asu-tsv) HTTP 接口，不需要注册任何账号。本模块把同一套接口封装成
MHP 认识的星表对象，并缓存同一视场的查询结果。

不要直接运行本文件；用法见 sip_calibrate.py --catalog。
"""

from __future__ import annotations

import re
import sys
import time
import urllib.parse
import urllib.request
from types import SimpleNamespace

import numpy as np

SOURCES = {
    "ucac4": {
        "vizier": "I/322A",
        "columns": "RAJ2000,DEJ2000,f.mag",
    },
    "gaia3": {
        "vizier": "I/355/gaiadr3",
        "columns": "RA_ICRS,DE_ICRS,Gmag",
    },
}

DEFAULT_MIRRORS = (
    "https://vizier.cds.unistra.fr",
    "https://vizier.u-strasbg.fr",
    "https://vizier.cfa.harvard.edu",
)

_ROW_RE = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([^\s]*)")


class VizieRCatalog:
    """Mimics the small Grappa3ECatalog interface used by MHP astrometry."""

    def __init__(
        self,
        catalog_root=None,
        source: str = "ucac4",
        mirrors=None,
        timeout: float = 120.0,
        max_rows: int = 200000,
        user_agent: str | None = None,
    ):
        if source not in SOURCES:
            raise ValueError(f"Unknown VizieR source: {source} (choose from {sorted(SOURCES)})")
        self.source = source
        self.mirrors = tuple(mirrors) if mirrors else DEFAULT_MIRRORS
        self.timeout = float(timeout)
        self.max_rows = int(max_rows)
        self.user_agent = (
            user_agent
            or "sip-align/1.0 (SIP calibration tool; VizieR public ASU)"
        )
        self._cache: dict[tuple, list[tuple[float, float, float]]] = {}

    def query_wcs_footprint(self, wcs, shape, margin_px: float = 45.0):
        """Download all catalog stars overlapping the WCS image footprint."""
        h, w = int(shape[0]), int(shape[1])
        m = float(margin_px)
        xs = np.array([-m, w - 1 + m, -m, w - 1 + m, (w - 1) / 2.0], dtype=float)
        ys = np.array([-m, -m, h - 1 + m, h - 1 + m, (h - 1) / 2.0], dtype=float)
        sky = wcs.pixel_to_world(xs, ys)
        center = sky[-1]
        sep_deg = float(np.max(center.separation(sky).deg))
        radius_arcmin = max(4.0, float(np.ceil(sep_deg * 60.0 + 3.0)))

        rows = self._fetch(center.ra.deg, center.dec.deg, radius_arcmin)
        if len(rows) < 15:
            raise RuntimeError(
                f"VizieR returned only {len(rows)} stars around "
                f"RA={center.ra.deg:.4f} Dec={center.dec.deg:+.4f} r={radius_arcmin:g}'"
            )

        ra = np.asarray([r[0] for r in rows], dtype=float)
        dec = np.asarray([r[1] for r in rows], dtype=float)
        mag = np.asarray([r[2] for r in rows], dtype=float)
        return SimpleNamespace(
            count=len(rows),
            files_read=[f"vizier:{self.source}"],
            ra_deg=ra,
            dec_deg=dec,
            g_mag_est=mag,
            # MHP only writes raw_value to its CSV as an integer.
            raw_value=np.zeros(len(rows), dtype=np.uint16),
        )

    def _fetch(self, ra_deg: float, dec_deg: float, radius_arcmin: float):
        key = (round(ra_deg, 5), round(dec_deg, 5), round(radius_arcmin, 2))
        cached = self._cache.get(key)
        if cached is not None:
            print(
                f"[VIZIER] cache hit {len(cached)} stars "
                f"r={radius_arcmin:g}' (same field)",
                file=sys.stderr,
            )
            return cached

        src = SOURCES[self.source]
        params = urllib.parse.urlencode(
            [
                ("-source", src["vizier"]),
                ("-c", f"{ra_deg:.6f} {dec_deg:+.6f}"),
                ("-c.rm", f"{radius_arcmin:g}"),
                ("-out.max", str(self.max_rows)),
                ("-out", src["columns"]),
            ]
        )
        last_exc: Exception | None = None
        for i, mirror in enumerate(self.mirrors):
            url = f"{mirror}/viz-bin/asu-tsv?{params}"
            try:
                req = urllib.request.Request(url, headers={"User-Agent": self.user_agent})
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    text = resp.read().decode("utf-8", errors="replace")
                rows = self._parse(text)
                if rows:
                    self._cache[key] = rows
                    return rows
                last_exc = RuntimeError(
                    f"{mirror}: 0 rows for RA={ra_deg:.4f} Dec={dec_deg:+.4f} r={radius_arcmin:g}'"
                )
            except Exception as exc:  # noqa: BLE001 - try next mirror
                last_exc = exc
                print(
                    f"[VIZIER] {mirror} failed: {type(exc).__name__}: {exc}",
                    file=sys.stderr,
                )
            if i + 1 < len(self.mirrors):
                time.sleep(1.0)
        raise RuntimeError(f"All VizieR mirrors failed: {last_exc}")

    @staticmethod
    def _parse(text: str):
        rows = []
        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = _ROW_RE.match(line)
            if not m:
                continue
            try:
                ra = float(m.group(1))
                dec = float(m.group(2))
            except ValueError:
                continue
            try:
                mag = float(m.group(3))
            except ValueError:
                mag = np.nan
            rows.append((ra, dec, mag))
        return rows
