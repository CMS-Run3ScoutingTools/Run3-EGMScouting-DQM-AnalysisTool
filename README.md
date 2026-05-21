# Run3-EGMScouting-DQM-AnalysisTool

Config-driven DQMIO pipeline with modular workflows:
- `mass_fit`: Crystal-Ball + background fits for invariant-mass histograms
- `tnp_efficiency`: TnP efficiency plots (era overlays, optional MC)
- `quantity_per_era`: per-era overlays of TnP quantity shapes using ROOT/cmsstyle
- `dqm_monitoring`: DQM `ScoutingMonitoring/Run summary` monitoring overlays and PAT/scouting comparisons
- `filter_monitoring`: filter-by-filter TnP monitoring using ROOT/cmsstyle canvases

## Architecture

- Core runner: `scripts/run_dqm_pipeline.py`
- Shared pipeline core (input resolution, DAS/golden/brilcalc, histogram loading/aggregation): `scripts/dqm_pipeline/core.py`
- Monitoring modules:
  - `scripts/dqm_pipeline/modules/mass_fit.py`
  - `scripts/dqm_pipeline/modules/tnp_efficiency.py`
  - `scripts/dqm_pipeline/modules/dqm_monitoring.py`
  - `scripts/dqm_pipeline/modules/filter_monitoring.py`
  - `scripts/dqm_pipeline/modules/quantity_per_era.py`
- Module registry: `scripts/dqm_pipeline/modules/__init__.py`

To add future work, implement a new module file with:

```python
def run_module(cfg, era_sources, out_root, strict=False):
    ...
    return summary_dict
```

Then register it in `scripts/dqm_pipeline/modules/__init__.py`.

## Environment

Load project environment (LCG + PYTHONPATH):

```bash
source env.sh
```

`env.sh` sources:

```bash
/cvmfs/sft.cern.ch/lcg/views/LCG_105b/x86_64-el9-gcc13-opt/setup.sh
```

After loading, run pipeline commands normally.

## Input model (new)

`eras` are resolved from DAS datasets (no hardcoded `file`, no fixed `run_number`):

- `DAS`: dataset path, e.g. `/ScoutingPFMonitor/Run2025G-PromptReco-v1/DQMIO`
- `file_glob`: local filesystem wildcard for EOS/DQMGUI copies, e.g. `/eos/.../DQM_*__Run2026B-PromptReco-v1__DQMIO.root`
- `run-requirement` (optional): supports `before`, `after`, `min`, `max`, `include`, `exclude`
- `golden_json`: local JSON path or `http(s)` URL for run filtering, configured per era
- `label` / `year` (optional): display metadata for mixed-year comparisons such as `Run2025G`, `Run2026A`, `Run2026B`

The script extracts run numbers from DQMIO filenames, filters runs, and aggregates histograms across selected runs per era.

For RelVal-style DQMIO files whose filenames do not encode the run number, enable preprocessing and set `run_number` explicitly. The preprocessing step converts DQMIO into a legacy DQM ROOT file that the rest of the framework can read normally.

Example preprocessing fields:

```yaml
preprocess_dqmio:
  enabled: true
  output_subdir: preprocessed_dqmio
  cmsrun_bin: cmsRun
  convention: Offline
  workflow_template: /RelVal/ScoutingPFMonitor/{era}
  skip_existing: true

eras:
  baseline:
    input_format: dqmio
    run_number: 333334
    file_glob: /eos/cms/store/relval/.../*.root
```

For cross-year comparisons, put all years in the same `eras:` block and assign each one its own `label`, `year`, and `golden_json`.
See `config/dqm_pipeline_mixed_years.example.yaml`.

## Golden JSON and lumi

Golden JSON is used for run filtering directly, and is configured per era.

If `lumi.use_brilcalc: true`, `lumi_fb` is estimated automatically per era using:
- selected runs from DAS + run filters
- golden JSON
- `brilcalc` output (`--output-style csv`)
- clean subprocess env by default (`lumi.clean_env: true`) to avoid Python env conflicts

When an era `golden_json` is a URL, the pipeline downloads it and materializes a local temp JSON automatically for `brilcalc -i`.

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

Example mixed-year plotting fields:

```yaml
plotting:
  campaign_label: Run 3
  lumi_text: Run 3 (13.6 TeV)
  energy_tev: 13.6
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

Run preprocessing only:

```bash
python3 scripts/run_dqm_pipeline.py --config config/dqm_comparison.yaml --preprocess-only
python3 scripts/preprocess_dqmio.py --config config/dqm_comparison.yaml
```

If `rich` is installed, running will show progress bars automatically.
Disable it with:

```bash
python3 scripts/run_dqm_pipeline.py --config config/dqm_pipeline.yaml --no-progress
```

Run one module:

```bash
python3 scripts/run_dqm_pipeline.py --config config/dqm_pipeline.yaml --module mass_fit
python3 scripts/run_dqm_pipeline.py --config config/dqm_pipeline.yaml --module tnp_efficiency
python3 scripts/run_dqm_pipeline.py --config config/dqm_pipeline.yaml --module dqm_monitoring
python3 scripts/run_dqm_pipeline.py --config config/dqm_pipeline.yaml --module quantity_per_era
python3 scripts/run_dqm_pipeline.py --config config/dqm_pipeline.yaml --module filter_monitoring
```

`filter_monitoring` supports:
- per-era filter-chain monitoring canvases
- era-by-era comparison for each filter (`save_era_comparison: true`)
- optional run-ordered auxiliary trend plots (`save_run_trends: true`)

Strict mode (fail immediately on missing file/hist):

```bash
python3 scripts/run_dqm_pipeline.py --config config/dqm_pipeline.yaml --strict
```
