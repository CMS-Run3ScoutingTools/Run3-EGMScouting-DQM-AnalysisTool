from dqm_pipeline.modules.mass_fit import run_module as run_mass_fit
from dqm_pipeline.modules.tnp_efficiency import run_module as run_tnp_efficiency


MODULE_REGISTRY = {
    "mass_fit": run_mass_fit,
    "tnp_efficiency": run_tnp_efficiency,
}
