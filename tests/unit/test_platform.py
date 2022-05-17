import os
import sys

from psyplot.config.logsetup import _get_home


def test_psy_stuff(name='psyplot', env_key='PSYPLOTCONFIGDIR'):
    print("Platform is dumdum...", sys.platform)
    p = None
    h = _get_home()
    print("Home is dumdum...", h)
    if ((sys.platform.startswith('linux') or sys.platform == 'darwin') and
            h is not None):
        p = os.path.join(h, '.config/' + name)
    elif h is not None:
        p = os.path.join(h, '.' + name)

    if not os.path.exists(p):
        os.makedirs(p)
