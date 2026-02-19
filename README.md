# Run3-EGMScouting-DQM-AnalysisTool

Config-driven DQMIO pipeline with modular workflows:
- `mass_fit`: Crystal-Ball + background fits for invariant-mass histograms
- `tnp_efficiency`: TnP efficiency plots (era overlays, optional MC)

## Architecture

- Core runner: `scripts/run_dqm_pipeline.py`
- Shared pipeline core (input resolution, DAS/golden/brilcalc, histogram loading/aggregation): `scripts/dqm_pipeline/core.py`
- Monitoring modules:
  - `scripts/dqm_pipeline/modules/mass_fit.py`
  - `scripts/dqm_pipeline/modules/tnp_efficiency.py`
- Module registry: `scripts/dqm_pipeline/modules/__init__.py`

To add future work, implement a new module file with:

```python
def run_module(cfg, era_sources, out_root, strict=False):
    ...
    return summary_dict
```

Then register it in `scripts/dqm_pipeline/modules/__init__.py`.

## Input model (new)

`eras` are resolved from DAS datasets (no hardcoded `file`, no fixed `run_number`):

- `DAS`: dataset path, e.g. `/ScoutingPFMonitor/Run2025G-PromptReco-v1/DQMIO`
- `run-requirement` (optional): supports `before`, `after`, `min`, `max`, `include`, `exclude`
- `golden_json`: local JSON for run filtering (global or era-specific)

The script extracts run numbers from DQMIO filenames, filters runs, and aggregates histograms across selected runs per era.

## Golden JSON and lumi

Golden JSON is used for run filtering directly.

If `lumi.use_brilcalc: true`, `lumi_fb` is estimated automatically per era using:
- selected runs from DAS + run filters
- golden JSON
- `brilcalc` output (`--output-style csv`)

## Mass Window Binning

In `mass_fit.fit_windows.<window>` you can now set either:
- `rebin`: integer rebin factor
- `nbins`: target number of bins inside that mass window

If both are set, `nbins` is used.

`mass_fit.overlay_rebin` controls only the era overlay plot binning.

Example config fields:

```yaml
lumi:
  use_brilcalc: true
  brilcalc_env: /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env
  brilcalc_bin: brilcalc
  calibration: web
  unit: /pb
```

## Run

Resolve eras/runs/files only:

```bash
python3 scripts/run_dqm_pipeline.py --config config/dqm_pipeline.yaml --resolve-only
```

Run all enabled modules:

```bash
python3 scripts/run_dqm_pipeline.py --config config/dqm_pipeline.yaml
```

Run one module:

```bash
python3 scripts/run_dqm_pipeline.py --config config/dqm_pipeline.yaml --module mass_fit
python3 scripts/run_dqm_pipeline.py --config config/dqm_pipeline.yaml --module tnp_efficiency
```

Strict mode (fail immediately on missing file/hist):

```bash
python3 scripts/run_dqm_pipeline.py --config config/dqm_pipeline.yaml --strict
```
