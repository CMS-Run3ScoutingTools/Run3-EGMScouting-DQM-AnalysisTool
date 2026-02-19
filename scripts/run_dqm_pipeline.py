#!/usr/bin/env python3
import argparse
from contextlib import contextmanager
from pathlib import Path

import ROOT
import yaml
import mplhep as hep

from dqm_pipeline.core import build_output_root, prepare_era_sources
from dqm_pipeline.modules import MODULE_REGISTRY

try:
    from rich.progress import (
        BarColumn,
        MofNCompleteColumn,
        Progress,
        TaskProgressColumn,
        TextColumn,
        TimeElapsedColumn,
        TimeRemainingColumn,
    )
except Exception:
    Progress = None


def parse_args():
    parser = argparse.ArgumentParser(description="Run modular DQMIO analysis pipeline from YAML config.")
    parser.add_argument("-c", "--config", required=True, help="Path to YAML config.")
    parser.add_argument(
        "--module",
        action="append",
        choices=[*MODULE_REGISTRY.keys(), "all"],
        help="Run only selected module(s). Can be repeated.",
    )
    parser.add_argument("--strict", action="store_true", help="Fail immediately on missing input/hist.")
    parser.add_argument("--resolve-only", action="store_true", help="Resolve eras/runs/files and exit.")
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
    if "resonance" not in cfg:
        raise RuntimeError("Config requires 'resonance'.")

    return cfg


def resolve_modules(requested):
    modules = set(requested or MODULE_REGISTRY.keys())
    if "all" in modules:
        modules = set(MODULE_REGISTRY.keys())
    return sorted(modules)


@contextmanager
def progress_context(enabled=True):
    use_rich = enabled and (Progress is not None)
    if not use_rich:
        yield None
        return

    with Progress(
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

    hep.style.use("CMS")

    out_root = build_output_root(cfg)
    out_root.mkdir(parents=True, exist_ok=True)

    with progress_context(enabled=not args.no_progress) as progress:
        era_sources = prepare_era_sources(
            cfg=cfg,
            config_dir=config_dir,
            strict=args.strict,
            progress=progress,
        )

        if args.resolve_only:
            print("[pipeline] resolve-only mode done.")
            return

        combined_summary = {}
        modules_to_run = resolve_modules(args.module)

        module_task = None
        if progress is not None:
            module_task = progress.add_task("[green]Running modules", total=len(modules_to_run))

        for module_name in modules_to_run:
            runner = MODULE_REGISTRY[module_name]
            combined_summary[module_name] = runner(
                cfg=cfg,
                era_sources=era_sources,
                out_root=out_root,
                strict=args.strict,
                progress=progress,
            )
            if progress is not None and module_task is not None:
                progress.update(module_task, advance=1)

    summary_file = out_root / "pipeline_summary.yaml"
    with open(summary_file, "w", encoding="utf-8") as f:
        yaml.safe_dump(combined_summary, f, sort_keys=False)

    print(f"[pipeline] Done. Main output directory: {out_root}")


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    main()
