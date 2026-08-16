"""Minimal path helper for the sip-align tool.

sip-align is usually run from a workspace that also contains the MHP project
(sibling directory). This module keeps the ``sys.path`` setup in one place so
the scripts can import ``m31_hmtproject`` helpers without hardcoded paths.
"""

from __future__ import annotations

import sys
from pathlib import Path


def ensure_on_sys_path(*roots: Path) -> None:
    """Insert project roots into ``sys.path`` once, at the front."""
    for root in roots:
        text = str(root)
        if text not in sys.path:
            sys.path.insert(0, text)


if __name__ == "__main__":
    print("sip-align path helper (ensure_on_sys_path)")
