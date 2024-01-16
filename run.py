from esmvalcore._main import run
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", module="iris")
    run()