"""Cross-engine regression: the full-LLG vgsot_llg.va must reproduce the Python
switching_vector trajectory. Skips gracefully when OpenVAF/ngspice are absent.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "va" / "llg" / "cross_validate.py"
SUMMARY = REPO / "va" / "llg" / "cross_validate_summary.json"


def _have(env_key, names):
    v = os.environ.get(env_key)
    if v and Path(v).exists():
        return True
    return any(shutil.which(n) for n in names)


_TOOLS = _have("OPENVAF", ["openvaf", "openvaf-r"]) and _have("NGSPICE", ["ngspice_con", "ngspice"])


@pytest.mark.skipif(not _TOOLS, reason="OpenVAF/ngspice not installed")
def test_vgsot_llg_va_matches_python_engine():
    r = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + "\n" + r.stderr
    s = json.loads(SUMMARY.read_text())
    assert s["pass"] is True
    assert s["max_abs_dmz"] < s["tol"]
    assert s["equilibrium_abs_dmz"] < 1e-3       # exact agreement at equilibrium
