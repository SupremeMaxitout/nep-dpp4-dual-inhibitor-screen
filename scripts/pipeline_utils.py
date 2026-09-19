#!/usr/bin/env python3
"""
Shared helpers used by every pipeline stage.

Nothing here is scientific. It is the plumbing: locating the repository root,
reading config/targets.yaml, setting up logging, finding external binaries, and
writing provenance records alongside every output.

The provenance part matters more than it looks. Every table this pipeline writes
gets a sibling .meta.json recording the config that produced it, the git commit,
the timestamp, and the versions of the tools involved. Six months from now that
is the difference between a result you can defend and a result you can only
apologise for. This repository exists because that record was not kept the first
time around.
"""

from __future__ import annotations

import json
import logging
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = REPO_ROOT / "config" / "targets.yaml"

DATA_RAW = REPO_ROOT / "data" / "raw"
DATA_INTERIM = REPO_ROOT / "data" / "interim"
DATA_PROCESSED = REPO_ROOT / "data" / "processed"
RESULTS = REPO_ROOT / "results"

for _d in (DATA_RAW, DATA_INTERIM, DATA_PROCESSED, RESULTS):
    _d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

class ConfigError(RuntimeError):
    """Raised when config/targets.yaml is missing something a stage needs."""


def load_config(path: Path | None = None) -> dict:
    """Read config/targets.yaml and refuse to proceed if placeholders remain."""
    path = path or CONFIG_PATH
    if not path.exists():
        raise ConfigError(f"Config not found at {path}")
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    if not cfg:
        raise ConfigError(f"Config at {path} is empty")
    return cfg


def get_target(cfg: dict, key: str) -> dict:
    """Fetch one target block, searching both targets and antitargets."""
    for section in ("targets", "antitargets"):
        if key in cfg.get(section, {}):
            block = dict(cfg[section][key])
            block["_key"] = key
            block["_section"] = section
            return block
    known = list(cfg.get("targets", {})) + list(cfg.get("antitargets", {}))
    raise ConfigError(f"Unknown target '{key}'. Known targets: {', '.join(known)}")


def check_no_placeholders(block: dict, key: str) -> None:
    """
    Walk a config block and fail loudly on any remaining VERIFY placeholder.

    This is deliberately fatal rather than a warning. A pipeline that runs to
    completion against an unverified structure produces results that look
    perfectly normal and are worthless.
    """
    bad: list[str] = []

    def walk(node: Any, trail: str) -> None:
        if isinstance(node, dict):
            for k, v in node.items():
                walk(v, f"{trail}.{k}" if trail else k)
        elif isinstance(node, list):
            for i, v in enumerate(node):
                walk(v, f"{trail}[{i}]")
        elif isinstance(node, str) and "VERIFY" in node.upper():
            bad.append(trail)

    walk(block, key)
    if bad:
        raise ConfigError(
            "Unresolved placeholders in config/targets.yaml:\n  "
            + "\n  ".join(bad)
            + "\n\nThese must be confirmed against the RCSB PDB before this stage "
              "can run. See docs/target_selection.md."
        )


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logging(name: str, verbose: bool = False) -> logging.Logger:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-7s  %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stdout,
    )
    return logging.getLogger(name)


# ---------------------------------------------------------------------------
# External tools
# ---------------------------------------------------------------------------

def require_binary(name: str, hint: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(
            f"Required binary '{name}' was not found on PATH.\n  {hint}\n"
            f"  Run  python scripts/00_check_environment.py  for a full report."
        )
    return path


def tool_version(binary: str) -> str:
    """Best-effort version string, for the provenance record."""
    if shutil.which(binary) is None:
        return "not installed"
    for flag in ("--version", "-v", "--help"):
        try:
            out = subprocess.run(
                [binary, flag], capture_output=True, text=True, timeout=10
            )
            text = (out.stdout or out.stderr).strip()
            if text:
                return text.splitlines()[0][:120]
        except Exception:
            continue
    return "unknown"


def git_commit() -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=REPO_ROOT, capture_output=True, text=True, timeout=10,
        )
        return out.stdout.strip() or "not a git repository"
    except Exception:
        return "unavailable"


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def write_provenance(output_path: Path, stage: str, params: dict,
                     tools: list[str] | None = None) -> Path:
    """
    Write <output>.meta.json next to any artefact the pipeline produces.

    Call this every time a stage writes a result. It costs nothing and it is
    the only reason anyone, including you, will be able to trust the output
    later.
    """
    meta = {
        "stage": stage,
        "output": output_path.name,
        "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": git_commit(),
        "python": sys.version.split()[0],
        "parameters": params,
        "tool_versions": {t: tool_version(t) for t in (tools or [])},
    }
    meta_path = output_path.with_suffix(output_path.suffix + ".meta.json")
    meta_path.write_text(json.dumps(meta, indent=2, default=str))
    return meta_path


def banner(logger: logging.Logger, title: str) -> None:
    logger.info("=" * 70)
    logger.info(title)
    logger.info("=" * 70)
