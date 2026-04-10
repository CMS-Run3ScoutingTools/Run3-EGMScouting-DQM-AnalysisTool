#!/usr/bin/env python3
import argparse
from contextlib import contextmanager
from pathlib import Path

import yaml

from dqm_pipeline.core import build_output_root, emit_log
from dqm_pipeline.preprocess_dqmio import preprocess_dqmio_inputs

try:
    from rich.progress import (
        BarColumn,
        MofNCompleteColumn,
        Progress,
        SpinnerColumn,
        TaskProgressColumn,
        TextColumn,
        TimeElapsedColumn,
        TimeRemainingColumn,
    )
except Exception:
    Progress = None


def parse_args():
    parser = argparse.ArgumentParser(description="Preprocess DQMIO inputs into legacy DQM ROOT files from YAML config.")
    parser.add_argument("-c", "--config", required=True, help="Path to YAML config.")
    parser.add_argument("--strict", action="store_true", help="Fail immediately on missing input/conversion failure.")
    parser.add_argument("--no-progress", action="store_true", help="Disable rich progress bar.")
    return parser.parse_args()


def load_config(config_path):
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise RuntimeError("Config must be a YAML mapping.")
    if "eras" not in cfg or not cfg["eras"]:
        raise RuntimeError("Config requires non-empty 'eras' section.")
    if "output_dir" not in cfg:
        raise RuntimeError("Config requires 'output_dir'.")
    return cfg


@contextmanager
def progress_context(enabled=True):
    use_rich = enabled and (Progress is not None)
    if not use_rich:
        if enabled and Progress is None:
            print("[preprocess][WARN] rich is not available; using plain text progress output.")
        yield None
        return

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
    ) as progress:
        yield progress


def main():
    args = parse_args()
    config_path = Path(args.config).resolve()
    config_dir = config_path.parent
    cfg = load_config(config_path)
    out_root = build_output_root(cfg)
    out_root.mkdir(parents=True, exist_ok=True)

    print(f"[preprocess] config={config_path}")
    print(f"[preprocess] out_root={out_root}")

    with progress_context(enabled=not args.no_progress) as progress:
        emit_log(progress, "[preprocess] stage: preprocess DQMIO", style="bold cyan")
        preprocess_dqmio_inputs(
            cfg=cfg,
            config_dir=config_dir,
            strict=args.strict,
            progress=progress,
        )
        emit_log(progress, "[preprocess] done.", style="bold green")


if __name__ == "__main__":
    main()
