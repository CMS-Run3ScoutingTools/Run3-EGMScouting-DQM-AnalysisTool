from dqm_pipeline.modules.dqm_monitoring import run_module as run_dqm_monitoring
from dqm_pipeline.modules.filter_monitoring import run_module as run_filter_monitoring
from dqm_pipeline.modules.mass_fit import run_module as run_mass_fit
from dqm_pipeline.modules.quantity_per_era import run_module as run_quantity_per_era
from dqm_pipeline.modules.tnp_efficiency import run_module as run_tnp_efficiency


MODULE_REGISTRY = {
    "dqm_monitoring": run_dqm_monitoring,
    "filter_monitoring": run_filter_monitoring,
    "mass_fit": run_mass_fit,
    "quantity_per_era": run_quantity_per_era,
    "tnp_efficiency": run_tnp_efficiency,
}
