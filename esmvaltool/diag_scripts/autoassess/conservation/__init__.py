"""Import conservation."""
from . import atmospheric_mass
from . import energy_budget
from . import GC_water_conservation

metrics_functions = [
    atmospheric_mass.global_atmos_mass_conservation,
    energy_budget.atmos_energy_budget,
    GC_water_conservation.global_freshwater_fluxes_over_various_GC_cubmodels
]
