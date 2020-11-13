"""Import stratosphere module and its functions."""
from . import age_of_air
from . import strat_metrics_1

metrics_functions = [strat_metrics_1.mainfunc, age_of_air.age_of_air]

multi_functions = [strat_metrics_1.multi_qbo_plot,
                   strat_metrics_1.multi_teq_plot,
                   strat_metrics_1.multi_t100_vs_q70_plot,
                   age_of_air.multi_age_plot]
